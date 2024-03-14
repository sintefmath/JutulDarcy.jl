function Jutul.prepare_step_storage(p::PrepareStepWellSolver, storage, model::MultiModel)
    @assert !isnothing(model.groups)
    submodels = Jutul.submodel_symbols(model)
    if false
        for (k, m) in pairs(model.models)
            if model_or_domain_is_well(m)
                ctrl_name = Symbol("$(k)_ctrl")
                ix = findfirst(isequal(ctrl_name), submodels)
                if isnothing(ix)
                    ix = findfirst(isequal(:Facility), submodels)
                end
                @assert !isnothing(ix) "No wells for well solver?"
                facility_label = submodels[ix]
                facility = model.models[facility_label]
                pos = findfirst(isequal(k), facility.domain.well_symbols)
            end
        end
    end

    targets = setdiff(submodels, (:Reservoir, ))
    groups = Int[]
    for (g, m) in zip(model.groups, submodels)
        if m in targets
            push!(groups, g)
        end
    end
    unique!(groups)
    well_solver_storage = JutulStorage()
    primary = JutulStorage()
    for target in targets
        primary[target] = deepcopy(storage[target].primary_variables)
    end
    well_solver_storage[:targets] = targets
    well_solver_storage[:well_primary_variables] = primary
    well_solver_storage[:groups_and_linear_solvers] = [(g, LUSolver()) for g in groups]
    return well_solver_storage
end

function Jutul.simulator_config!(cfg, sim, ::PrepareStepWellSolver, storage)
    add_option!(cfg, :well_iterations, 25, "Well iterations to be performed before step.", types = Int, values = 0:10000)
    add_option!(cfg, :well_outer_iterations_limit, Inf, "Number of outer iterations for each solve where wells are to be solved.", types = Union{Int, Float64})
    add_option!(cfg, :well_acceptance_factor, 10.0, "Accept well pre-solve results at this relaxed factor.", types = Float64)
    add_option!(cfg, :well_info_level, -1, "Info level for well solver.", types = Int)
end

function Jutul.prepare_step!(wsol_storage, wsol::PrepareStepWellSolver, storage, model::MultiModel, dt, forces, config;
        executor = DefaultExecutor(),
        iteration = 0,
        relaxation = 1.0
    )
    if iteration <= config[:well_outer_iterations_limit]
        converged = false
        num_its = -1
        max_well_iterations = config[:well_iterations]
        targets = wsol_storage.targets
        il = config[:well_info_level]

        primary = wsol_storage.well_primary_variables
        for target in targets
            for (k, v) in pairs(storage[target].primary_variables)
                primary[target][k] .= v
            end
        end
        converged, num_its = solve_well_system_with_fixed_reservoir(
            storage, model, dt, forces, executor,
            wsol_storage.groups_and_linear_solvers, targets, iteration, relaxation,
            il, max_well_iterations, config)
        if converged
            if il > 0
                jutul_message("Well solver", "Converged in $num_its iterations.")
            end
        else
            # TODO: Should only cover non-converged wells.
            for target in targets
                for (k, v) in pairs(storage[target].primary_variables)
                    v .= primary[target][k]
                end
            end
            Jutul.update_state_dependents!(storage, model, dt, forces,
                update_secondary = true, targets = targets)
            if il > 0
                jutul_message("Well solver", "Did not converge in $max_well_iterations iterations.", color = :yellow)
            end
        end
    end

    return (nothing, forces)
end

function solve_well_system_with_fixed_reservoir(storage, model, dt, forces, executor, groups_and_lsolve, targets, iteration, relaxation, il, max_well_iterations, config)
    for well_it in 1:max_well_iterations
        Jutul.update_state_dependents!(storage, model, dt, forces,
            update_secondary = well_it > 1,
            targets = targets
        )
        Jutul.update_linearized_system!(storage, model, executor, targets = targets)
        if well_it == max_well_iterations
            tol_factor = config[:well_acceptance_factor]
        else
            tol_factor = 1.0
        end
        ok = map(
            k -> check_convergence(storage[k], model[k], config[:tolerances][k],
                iteration = iteration,
                dt = dt,
                tol_factor = tol_factor
            ),
            targets
        )
        if il > 0
            jutul_message("Well solver #$well_it/$max_well_iterations", "$(sum(ok))/$(length(ok)) well and facility models are converged.", color = :cyan)
            if il > 1
                for (k, ok_i) in zip(targets, ok)
                    if !ok_i
                        jutul_message("$k", "Not converged", color = :dark_gray)
                    end
                end
            end
        end
        if all(ok)
            return (true, well_it)
        end
        for (g, lsolve) in groups_and_lsolve
            lsys = storage.LinearizedSystem[g, g]
            recorder = storage.recorder
            linear_solve!(lsys, lsolve, model, storage, dt, recorder, executor)
        end
        update = Jutul.update_primary_variables!(storage, model; targets = targets, relaxation = relaxation)
    end
    return (false, max_well_iterations)
end
