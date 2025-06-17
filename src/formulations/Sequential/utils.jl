

function setup_sequential_storage!(S, p_model::MultiModel, t_model::MultiModel)
    for k in Jutul.submodels_symbols(p_model)
        p_submodel = p_model[k]
        t_submodel = t_model[k]
        if k == :Reservoir || model_or_domain_is_well(p_submodel)
            S[k] = JutulStorage()
            setup_sequential_storage!(S[k], p_submodel, t_submodel)
        end
    end
end

function setup_sequential_storage!(S, p_model::SimulationModel, t_model::SimulationModel)
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
    transfer_keys_transport = transfer_keys(t_model, p_model)
    transfer_keys_pressure = transfer_keys(p_model, t_model)
    # Init rule: Transfer over everything that is present from transport state.
    function init_keys(target, source)
        all_t = keys(Jutul.get_variables_by_type(target, :all))
        all_s = keys(Jutul.get_variables_by_type(source, :all))

        # return intersect(all_t, all_s)

        secondary_s = keys(Jutul.get_variables_by_type(source, :secondary))
        tkeys = setdiff(intersect(all_t, all_s), secondary_s)
        push!(tkeys, :TotalMasses)
        return tkeys
    end
    init_keys_pressure = init_keys(p_model, t_model)

    transfer_keys = JutulStorage()
    transfer_keys[:pressure] = transfer_keys_pressure
    transfer_keys[:transport] = transfer_keys_transport

    init_keys = JutulStorage()
    init_keys[:pressure] = init_keys_pressure
    # init_keys[:transport] = init_keys_transport

    S[:init_keys] = init_keys
    S[:transfer_keys] = transfer_keys

    nph = number_of_phases(t_model.system)
    nc = number_of_cells(t_model.domain)
    λ = zeros(nph, nc)
    @. λ = NaN
    S[:mobility] = λ
    S[:mobility_prev] = similar(λ)
    return S
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

function Jutul.simulator_executor(sim::SequentialSimulator)
    return Jutul.simulator_executor(sim.pressure)
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

function get_reservoir_state(sim, current = true)
    if current
        model_state = sim.storage.state
    else
        model_state = sim.storage.state0
    end
    if sim.model isa MultiModel
        model_state = model_state.Reservoir
    end
    return model_state
end

function Jutul.perform_step!(
        simulator::SequentialSimulator,
        dt,
        forces,
        config;
        iteration = NaN,
        relaxation = 1.0,
        update_secondary = true,
        solve = true,
        executor = default_executor(),
        prev_report = missing
    )
    il = config[:info_level]

    psim = simulator.pressure
    pstate = psim.storage.state
    pstate0 = psim.storage.state0
    storage = simulator.storage

    is_multi_model = psim.model isa MultiModel
    tsim = simulator.transport
    tstate = tsim.storage.state

    # Solve pressure
    # Copy over variables to parameters for both solves
    if iteration > 1
        # We need to transfer from pressure
        sequential_sync_values!(simulator, to_key = :pressure)
    end
    report = Jutul.setup_ministep_report()
    config_p = config[:pressure]
    max_iter_p = config_p[:max_nonlinear_iterations]
    mob_p = get_reservoir_state(psim).PhaseMobilities
    mob_t = get_reservoir_state(tsim).PhaseMobilities
    total_saturation = get_reservoir_state(tsim).TotalSaturation

    transfer_storage = simulator.storage.transfer
    if is_multi_model
        mob = transfer_storage.Reservoir.mobility
        mob_prev = transfer_storage.Reservoir.mobility_prev
    else
        mob = transfer_storage.mobility
        mob_prev = transfer_storage.mobility_prev
    end
    if iteration == 1
        if isnan(mob[1, 1])
            tstate0 = get_reservoir_state(tsim, false)
            mob_t0 = tstate0.PhaseMobilities
            @. mob = value(mob_t0)
        else
            # Initial guess from end of last time-step
            @. mob = value(mob_t)
        end
    end

    @. mob_p = mob
    @. mob_prev = mob
    if il > 1
        jutul_message("Sequential it $iteration", "Solving pressure")
    end
    done_p, report_p = Jutul.solve_ministep(psim, dt, forces, max_iter_p, config_p, finalize = false)
    if done_p
        # t_forces = transport_forces(tsim,)
        t_forces = deepcopy(forces)
        # Copy over values for pressure and fluxes into parameters for second simulator
        model_p = psim.model
        if iteration > 1
            Jutul.reset_previous_state!(tsim, pstate0)
            Jutul.reset_state_to_previous_state!(tsim)
        end
        store_fluxes!(tsim, psim)
        sequential_sync_values!(simulator, to_key = :transport)
        nsub = config[:transport_substeps]
        config_t = config[:transport]
        max_iter_t = config_t[:max_nonlinear_iterations]
        # TODO: Store initial guesses here for SFI

        function store_mobility!(mob, mob_t, w, cell_weight)
            cw(::Nothing, i) = 1.0
            cw(x, i) = value(x[i])
            for i in axes(mob_t, 2)
                # λ_t = 0.0
                # for ph in axes(mob_t, 1)
                #   λ_t += value(mob_t[ph, i])
                # end
                w_i = w*cw(cell_weight, i)
                for ph in axes(mob_t, 1)
                    λ = value(mob_t[ph, i])
                    # mob[ph, i] += λ/λ_t
                    mob[ph, i] += w_i*λ
                end
            end
        end
        report_t = nothing
        @. mob = 0
        if il > 1
            jutul_message("Sequential it $iteration", "Solving transport with $nsub substeps")
        end
        if nsub == 1
            # Then transport
            done_t, report_t = Jutul.solve_ministep(tsim, dt, t_forces, max_iter_t, config_t)
            store_mobility!(mob, mob_t, 1.0, nothing)
        else
            for stepno = 1:nsub
                dt_i = dt/nsub
                done_t, subreport_t = Jutul.solve_ministep(tsim, dt_i, t_forces, max_iter_t, config_t)
                if stepno == nsub
                    w_i = nothing
                else
                    w_i = total_saturation
                end
                store_mobility!(mob, mob_t, dt_i, w_i)
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
    else
        report[:failure] = true
        report_t = missing
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
            if ismissing(R)
                continue
            end
            for step in R[:steps]
                if haskey(step, k)
                    v += step[k]
                end
            end
        end
        report[k] = v
    end
    t_prep = 0.0
    t_precond = 0.0
    precond_count = 0
    for R in [report_p, report_t]
        if ismissing(R)
            continue
        end
        for step in R[:steps]
            lsol = get(step, :linear_solver, nothing)
            if !isnothing(lsol)
                t_prep += lsol.prepare
                t_precond += lsol.precond
                precond_count += lsol.precond_count
            end
        end
    end
    report[:linear_solver] = (
        stats = nothing,
        prepare = t_prep,
        precond = t_precond,
        precond_count = precond_count
    )
    # Return convergence criterion for outer loop if SFI
    sfi = config[:sfi]
    if !done_p
        report[:failure_exception] = ErrorException("Pressure solve failed to converge.")
        converged = false
        err = Inf
    elseif !done_t
        report[:failure_exception] = ErrorException("Transport solve failed to converge.")
        converged = false
        err = Inf
    elseif sfi
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
        for sT in total_saturation
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
        converged = true
        err = 0.0
    end
    report[:converged] = converged
    if converged
        Jutul.update_after_step!(psim, dt, forces)
    end
    return (err, converged, report)
end

function store_fluxes!(tsim, psim)
    pmodel = psim.model
    vT = get_reservoir_state(tsim).TotalVolumetricFlux
    pstate_res = get_reservoir_state(psim)
    store_total_fluxes!(vT, reservoir_model(pmodel), as_value(pstate_res))
    store_perforation_fluxes!(tsim, psim, pmodel)
end

function store_perforation_fluxes!(tsim, psim, pmodel)
    nothing
end

function store_perforation_fluxes!(tsim, psim, pmodel::MultiModel)
    model_r = reservoir_model(pmodel)
    for ct in pmodel.cross_terms
        if ct.cross_term isa PressureReservoirFromWellFlowCT
            @assert ct.target == :Reservoir "Expected target to be :Reservoir, was $(ct.target) (source = $(ct.source))"
            label = ct.source
            q = tsim.storage[label].parameters[:PerforationTotalVolumetricFlux]
            ct_mphase = ct.cross_term.parent
            store_perforation_fluxes!(q, ct_mphase, pmodel[label], model_r, psim.storage.state[label], psim.storage.state[:Reservoir])
        end
    end
    return tsim
end

function store_perforation_fluxes!(q, ct, wmodel, resmodel, wstate, rstate)
    T = eltype(q)
    nph = number_of_phases(wmodel.system)
    mob = rstate.PhaseMobilities
    for i in eachindex(q)
        conn = cross_term_perforation_get_conn(ct, i, wstate, rstate)
        total_mobility = zero(T)
        for ph in 1:nph
            total_mobility += value(mob[ph, conn.reservoir])
        end
        total_flux = zero(T)
        for ph in 1:nph
            pot = value(JutulDarcy.perforation_phase_potential_difference(conn, rstate, wstate, ph))
            total_flux += pot*total_mobility
        end
        q[i] = total_flux
    end
    return q
end

function sequential_sync_values!(sim::SequentialSimulator; to_key::Symbol = :pressure, kind = :transfer_keys)
    @assert to_key in (:pressure, :transport)
    pstate = sim.pressure.storage.state
    tstate = sim.transport.storage.state
    if to_key == :pressure
        dest = pstate
        src = tstate
        model = sim.transport.model
    else
        dest = tstate
        src = pstate
        model = sim.pressure.model
    end
    sequential_sync_values!(dest, src, model, sim.storage.transfer, to_key, kind)
end

function sequential_sync_values!(dest, src, model::SimulationModel, storage_transfer, dest_key, kind)
    for k in storage_transfer[kind][dest_key]
        update_values!(dest[k], src[k])
    end
end

function sequential_sync_values!(dest, src, model::MultiModel, storage_transfer, k, kind)
    for modname in keys(storage_transfer)
        sequential_sync_values!(dest[modname], src[modname], model[modname], storage_transfer[modname], k, kind)
    end
end

# function Jutul.update_after_step!(sim::SequentialSimulator, dt, forces; kwarg...)
#     report = Dict{Symbol, Any}()
#     report[:pressure] = Jutul.update_after_step!(sim.pressure, dt, forces; kwarg...)
#     report[:transport] = Jutul.update_after_step!(sim.transport, dt, forces; kwarg...)
#     return report
# end

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
                    if !haskey(step, k)
                        continue
                    end
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
    sequential_sync_values!(sim, to_key = :pressure, kind = :init_keys)
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
