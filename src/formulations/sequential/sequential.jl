struct SequentialSimulator{M, P, T, S} <: Jutul.JutulSimulator
    model::M
    pressure::P
    transport::T
    storage::S
end

function SequentialSimulator(case::JutulCase; kwarg...)
    return SequentialSimulator(case.model; state0 = case.state0, parameters = case.parameters)
end

function SequentialSimulator(model; state0 = setup_state(model), parameters = setup_parameters(model), avg_mobility = false)
    rmodel = reservoir_model(model)
    sys = rmodel.system
    pmodel = convert_to_sequential(model, pressure = true, avg_mobility = avg_mobility)
    tmodel = convert_to_sequential(model, pressure = false)
    function add_total_saturation!(m, state0)
        if !haskey(state0, :TotalSaturation)
            state0[:TotalSaturation] = ones(number_of_cells(m.domain))
        end
    end
    function add_total_saturation!(m::MultiModel, state0)
        add_total_saturation!(m[:Reservoir], state0[:Reservoir])
    end
    add_total_saturation!(model, state0)

    function merge_initial_state(m, state0, parameters)
        return merge(state0, parameters)
    end
    function merge_initial_state(m::MultiModel, state0, parameters)
        init = copy(state0)
        init[:Reservoir] = merge(state0[:Reservoir], parameters[:Reservoir])
        return init
    end
    init = merge_initial_state(model, state0, parameters)
    function subsimulator(m)
        s0, prm = setup_state_and_parameters(m, init)
        return Simulator(m, state0 = s0, parameters = prm)
    end

    PSim = subsimulator(pmodel)
    TSim = subsimulator(tmodel)

    function seq_output_keys(m)
        return copy(m.output_variables)
    end
    # Transfer rule: All parameters present in other that are not parameters both places.
    function transfer_keys(target, source)
        prm_t = keys(Jutul.get_parameters(target))
        prm_s = keys(Jutul.get_parameters(source))

        all_t = keys(Jutul.get_variables_by_type(target, :all))
        all_s = keys(Jutul.get_variables_by_type(source, :all))

        return intersect(all_s, setdiff(prm_t, prm_s))
    end
    transfer_keys_transport = transfer_keys(TSim.model, PSim.model)
    transfer_keys_pressure = transfer_keys(PSim.model, TSim.model)
    # Init rule: Transfer over everything that is present from transport state.
    function init_keys(target, source)
        all_t = keys(Jutul.get_variables_by_type(target, :all))
        all_s = keys(Jutul.get_variables_by_type(source, :all))

        return intersect(all_t, all_s)
    end
    # init_keys_transport = init_keys(TSim.model)
    init_keys_pressure = init_keys(PSim.model, TSim.model)

    transfer_keys = JutulStorage()
    transfer_keys[:pressure] = transfer_keys_pressure
    transfer_keys[:transport] = transfer_keys_transport

    init_keys = JutulStorage()
    init_keys[:pressure] = init_keys_pressure
    # init_keys[:transport] = init_keys_transport

    S = JutulStorage()
    S[:init_keys] = init_keys
    S[:transfer_keys] = transfer_keys

    nph = number_of_phases(sys)
    nc = number_of_cells(model.domain)
    λ = zeros(nph, nc)
    @. λ = NaN
    S[:mobility] = λ
    S[:mobility_prev] = similar(λ)

    # @info "Keys" transfer_keys_pressure transfer_keys_transport init_keys_transport init_keys_pressure


    # all_keys = union(seq_output_keys(PSim.model), seq_output_keys(TSim.model))
    # seq_state = JutulStorage()
    # for s in [PSim, TSim]
    #     s0 = s.storage.state0
    #     for k in seq_output_keys(s.model)
    #         if haskey(s0, k)
    #             seq_state[k] = similar(s0[k])
    #         end
    #     end
    # end
    # @info "!" seq_state
    # error()

    S[:recorder] = ProgressRecorder()
    # S[:state] = seq_state
    # S[:state0] = deepcopy(seq_state)
    # S[:primary_variables] = seq_state
    return SequentialSimulator(model, PSim, TSim, S)
end

function Jutul.simulator_config(
        sim::SequentialSimulator;
        transport_substeps = 1,
        saturation_tol = 1e-2,
        mobility_tol = 1e-2,
        sfi = true,
        kwarg...
    )
    cfg = Jutul.JutulConfig("Simulator config")
    Jutul.add_option!(cfg, :transport_substeps, transport_substeps, types = Int)
    Jutul.add_option!(cfg, :saturation_tol, saturation_tol, types = Float64)
    Jutul.add_option!(cfg, :mobility_tol, mobility_tol, types = Float64)
    Jutul.add_option!(cfg, :sfi, sfi, types = Bool)

    Jutul.simulator_config!(cfg, sim;
        min_nonlinear_iterations = 0,
        kwarg...
    )

    for k in [:pressure, :transport]
        cfg_k = Jutul.simulator_config(getfield(sim, k); always_update_secondary = true, kwarg...)
        Jutul.add_option!(cfg, k, cfg_k)
    end
    return cfg
end



function Jutul.select_linear_solver(sim::SequentialSimulator; kwarg...)
    return nothing
end

# function Jutul.reset_variables!(sim::SequentialSimulator, vars; kwarg...)
#     Jutul.reset_variables!(sim.pressure, vars; kwarg...)
#     Jutul.reset_variables!(sim.transport, vars; kwarg...)
#     return sim
# end

function Jutul.initial_setup!(sim::SequentialSimulator, config, timesteps; kwarg...)
    Jutul.initial_setup!(sim.pressure, config, timesteps; kwarg...)
    Jutul.initial_setup!(sim.transport, config, timesteps; kwarg...)
end

function Jutul.initialize_before_first_timestep!(sim::SequentialSimulator, dt; kwarg...)
    Jutul.initialize_before_first_timestep!(sim.pressure, dt; kwarg...)
end

function Jutul.perform_step!(
        simulator::SequentialSimulator,
        dt,
        forces,
        config;
        iteration = NaN,
        relaxation = 1.0,
        update_secondary = true,
        solve = true
    )
    psim = simulator.pressure
    pstate = psim.storage.state
    pstate0 = psim.storage.state0

    tsim = simulator.transport
    tstate = tsim.storage.state

    transfer_keys = simulator.storage.transfer_keys
    # Solve pressure
    # Copy over variables to parameters for both solves
    if iteration > 1
        # We need to transfer from pressure
        for k in transfer_keys[:pressure]
            update_values!(pstate[k], tstate[k])
        end
    end
    report = Jutul.setup_ministep_report()
    config_p = config[:pressure]
    max_iter_p = config_p[:max_nonlinear_iterations]
    mob_p = psim.storage.state.PhaseMobilities
    mob_t = tsim.storage.state.PhaseMobilities
    mob = simulator.storage.mobility
    mob_prev = simulator.storage.mobility_prev
    if iteration == 1
        if isnan(mob[1, 1])
            mob_t0 = tsim.storage.state0.PhaseMobilities
            @. mob = value(mob_t0)
        else
            # Initial guess from end of last time-step
            @. mob = value(mob_t)
        end
    end

    @. mob_p = mob
    @. mob_prev = mob

    done_p, report_p = Jutul.solve_ministep(psim, dt, forces, max_iter_p, config_p, finalize = false)
    if done_p
        # Copy over values for pressure and fluxes into parameters for second simulator
        model_p = psim.model
        if iteration > 1
            Jutul.reset_previous_state!(tsim, pstate0)
            Jutul.reset_state_to_previous_state!(tsim)
        end

        vT = tsim.storage.state.TotalVolumetricFlux
        store_total_fluxes!(vT, model_p, as_value(pstate))
        for k in transfer_keys[:transport]
            update_values!(tstate[k], pstate[k])
        end
        nsub = config[:transport_substeps]
        config_t = config[:transport]
        max_iter_t = config_t[:max_nonlinear_iterations]
        # TODO: Store initial guesses here for SFI

        function store_mobility!(mob, mob_t, w)
            for i in axes(mob_t, 2)
                # λ_t = 0.0
                # for ph in axes(mob_t, 1)
                #   λ_t += value(mob_t[ph, i])
                # end
                for ph in axes(mob_t, 1)
                    λ = value(mob_t[ph, i])
                    # mob[ph, i] += λ/λ_t
                    mob[ph, i] += w*λ
                end
            end
        end
        report_t = nothing
        @. mob = 0
        if nsub == 1
            # Then transport
            done_t, report_t = Jutul.solve_ministep(tsim, dt, forces, max_iter_t, config_t)
            store_mobility!(mob, mob_t, 1.0)
        else
            for stepno = 1:nsub
                dt_i = dt/nsub
                done_t, subreport_t = Jutul.solve_ministep(tsim, dt_i, forces, max_iter_t, config_t)
                store_mobility!(mob, mob_t, dt_i)
                if stepno > 1
                    for (k, v) in subreport_t
                        if k == :steps
                            for s in v
                                push!(report_t[k], s)
                            end
                        elseif v isa AbstractFloat
                            report_t[k] += v
                        else
                            # We skip some data.
                            report_t[k] = v
                        end
                    end
                else
                    report_t = subreport_t
                end
            end
            @. mob /= dt
        end
        report[:transport] = report_t

    else
        error("Pressure failure not implemented")
    end
    report[:pressure] = report_p
    for k in [
        :secondary_time,
        :equations_time,
        :linear_system_time,
        :convergence_time,
        :linear_solve_time,
        :update_time,
        :linear_iterations
        ]
        v = 0
        for R in [report_p, report_t]
            for step in R[:steps]
                if haskey(step, k)
                    v += step[k]
                end
            end
        end
        report[k] = v
    end
    # Return convergence criterion for outer loop if SFI
    sfi = config[:sfi]
    if sfi
        tol_s = 1e-2
        tol_mob = 1e-2
        tol_mob = config[:mobility_tol]
        tol_s = config[:saturation_tol]
        il = config[:info_level]
        min_its = config[:min_nonlinear_iterations]

        e_mob = 0
        for i in axes(mob, 2)
            λ_t = 0.0
            for ph in axes(mob_t, 1)
                λ_t += mob[ph, i]
            end
            for ph in axes(mob, 1)
                e = abs(mob[ph, i] - mob_prev[ph, i])/λ_t
                e_mob = max(e, e_mob)
            end
        end

        e_s = 0.0
        for sT in tsim.storage.state.TotalSaturation
            e_s = max(abs(value(sT) - 1.0), e_s)
        end

        converged = e_s <= tol_s && e_mob <= tol_mob
        converged = converged && iteration > min_its
        err = max(e_s, e_mob)

        # errors = [
        #     (
        #         name = :total_saturation,
        #         tolerances = Dict(:total_saturation => tol_s),
        #         criterions = ((
        #             E = (errors = e_s, ),
        #             names = (:E1, ),
        #         ),),
        #     )
        # ]
        if il > 1 # || true
            jutul_message("#$iteration", "|S_t - 1| = $e_s, |Δλ| = $e_mob")
           # Jutul.get_convergence_table(errors, il, iteration, config)
        end
    else
        converged = done_t
        err = 0.0
    end
    report[:converged] = converged
    if converged
        Jutul.update_after_step!(psim, dt, forces)
    end
    return (err, converged, report)
end

function Jutul.update_after_step!(sim::SequentialSimulator, dt, forces; kwarg...)
    report = Dict{Symbol, Any}()
    report[:pressure] = Jutul.update_after_step!(sim.pressure, dt, forces; kwarg...)
    report[:transport] = Jutul.update_after_step!(sim.transport, dt, forces; kwarg...)
    return report
end

function Jutul.get_output_state(sim::SequentialSimulator)
    Jutul.get_output_state(sim.transport)
end

function Jutul.final_simulation_message(simulator::SequentialSimulator, p, rec, t_elapsed, reports, arg...)
    Jutul.jutul_message("Sequential", "Total timing, per SFI iteration")
    Jutul.final_simulation_message(simulator.pressure, p, rec, t_elapsed, reports, arg...)

    # function make_sub_report(reports, typ)
    #     reports_seq = Vector{Any}()
    #     for rep in reports
    #         rep_seq = similar(rep)
    #         for (k, v) in rep
    #             if k == :ministeps
    #                 @info "?!" v
    #                 rep_inner = Vector{Any}()
    #                 for ministep in v
    #                     @info ministep rep_inner
    #                 end

    #                 rep_seq[k] = rep_inner
    #             else
    #                 rep_seq[k] = v
    #             end
    #         end
    #         push!(reports_seq, rep_seq)
    #     end
    #     return reports_seq
    # end
    # @info "Hey" make_sub_report(reports, :pressure)


    stats = Dict()
    stats[:pressure] = Dict(:iterations => 0, :linear_iterations => 0)
    stats[:transport] = Dict(:iterations => 0, :linear_iterations => 0)

    for rep in reports
        for ministep in rep[:ministeps]
            for step in ministep[:steps]
                for k in [:pressure, :transport]
                    for sstep in step[k][:steps]
                        if haskey(sstep, :update)
                            stats[k][:iterations] += 1
                            stats[k][:linear_iterations] += sstep[:linear_iterations]
                        end
                    end
                end
            end
        end
    end
    for k in [:pressure, :transport]
        s = stats[k]
        i = s[:iterations]
        li = s[:linear_iterations]
        Jutul.jutul_message(titlecase("$k"), "$i iterations, $li linear iterations")
    end
    # Jutul.jutul_message("Transport")
    # Jutul.final_simulation_message(simulator.transport, arg...)
end

function Jutul.update_before_step!(sim::SequentialSimulator, dt, forces; kwarg...)
    Jutul.update_before_step!(sim.pressure, dt, forces; kwarg...)
    Jutul.update_before_step!(sim.transport, dt, forces; kwarg...)
end

function Jutul.update_after_step!(sim::SequentialSimulator, dt, forces; kwarg...)
    # First, sync up state and state0 for transport
    Jutul.update_after_step!(sim.transport, dt, forces; kwarg...)
    # NOTE: This part must be done carefully so that the pressure contains the
    # final solution from the last transport solve.
    pstate = sim.pressure.storage.state
    tstate = sim.transport.storage.state
    for k in sim.storage.init_keys.pressure
        update_values!(pstate[k], tstate[k])
    end
    Jutul.update_after_step!(sim.pressure, dt, forces; kwarg...)
end


function Jutul.reset_state_to_previous_state!(sim::SequentialSimulator)
    Jutul.reset_state_to_previous_state!(sim.pressure)
    Jutul.reset_state_to_previous_state!(sim.transport)
end


function Jutul.reset_variables!(sim::SequentialSimulator, state; kwarg...)
    Jutul.reset_variables!(sim.pressure, state; kwarg...)
    Jutul.reset_variables!(sim.transport, state; kwarg...)
end

function Jutul.reset_previous_state!(sim::SequentialSimulator, state)
    Jutul.reset_previous_state!(sim.pressure, state)
    Jutul.reset_previous_state!(sim.transport, state)
end
# reset_previous_state!

include("interface.jl")
