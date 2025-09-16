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
        # t_forces = deepcopy(forces)
        # Copy over values for pressure and fluxes into parameters for second simulator
        # model_p = psim.model
        if iteration > 1
            Jutul.reset_previous_state!(tsim, pstate0)
            Jutul.reset_state_to_previous_state!(tsim)
        end
        store_fluxes!(tsim, psim)
        sequential_sync_values!(simulator, to_key = :transport)
        t_forces = transport_forces(tsim, psim, forces)
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
    report[:post_update_seq] = Dict{Symbol, Any}(
        :pressure => missing,
        :transport => missing
    )
    for (label, R) in [(:pressure, report_p), (:transport, report_t)]
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
        if haskey(R, :post_update)
            report[:post_update_seq][label] = R[:post_update]
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
        # Assume that simple wells are used (i.e. equal potential difference for all phases)
        WI_dp = JutulDarcy.perforation_phase_potential_difference(conn, rstate, wstate, 1)
        # Positive flux into reservoir
        q[i] = -value(WI_dp)*total_mobility
    end
    return q
end

function transport_forces(tsim, psim, forces)
    t_forces = deepcopy(forces)
    tmodel = tsim.model
    pmodel = psim.model
    fkey = :Facility

    has_facility = pmodel isa MultiModel && haskey(pmodel.models, :Facility)
    haskey(t_forces, fkey) == has_facility || error("Transport forces does not have key $fkey")

    if has_facility
        facility_model = tsim.model[fkey]
        fstate = psim.storage[fkey].state
        new_controls = Dict{Symbol, Any}()
        for (wno, wname) in enumerate(tsim.model[fkey].domain.well_symbols)
            # Note: To handle well switches we need to take the current control from the pressure solver.
            # ctrl = t_forces[fkey].control[wname]
            ctrl = fstate.WellGroupConfiguration.operating_controls[wname]
            wstate = psim.storage.state[wname]
            if ctrl isa JutulDarcy.InjectorControl
                if ctrl.target isa JutulDarcy.BottomHolePressureTarget
                    q_w = JutulDarcy.compute_well_qoi(pmodel, psim.storage.state, forces, wname, :rate)
                    q_w = Jutul.value(q_w)
                    new_ctrl = JutulDarcy.replace_target(ctrl, JutulDarcy.TotalRateTarget(q_w))
                else
                    new_ctrl = ctrl
                end
            elseif ctrl isa JutulDarcy.ProducerControl
                q_w_perf = compute_well_volumetric_rate_at_bhp(wstate, fstate.TotalSurfaceMassRate[wno])
                new_target = JutulDarcy.ReservoirVolumeRateTarget(q_w_perf)
                new_ctrl = JutulDarcy.replace_target(ctrl, new_target)
            elseif ctrl isa JutulDarcy.DisabledControl
                new_ctrl = ctrl
            else
                error("Unhandled control type $(typeof(ctrl)) for well $wname in transport_forces")
            end
            new_controls[wname] = new_ctrl
        end
        t_forces[fkey] = Jutul.setup_forces(facility_model, control = new_controls, set_default_limits = false)
    end
    return t_forces
end

function compute_well_volumetric_rate_at_bhp(wstate, well_mass_rate)
    rho = wstate.PhaseMassDensities
    s = wstate.Saturations
    total_density = 0.0
    c = JutulDarcy.well_top_node()
    for ph in axes(rho, 1)
        total_density += Jutul.value(rho[ph, c])*Jutul.value(s[ph, c])
    end
    q_vol = Jutul.value(well_mass_rate)/total_density
    q_vol::Float64
    return q_vol
end