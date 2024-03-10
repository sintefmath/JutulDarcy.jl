import Jutul: update_cross_terms!, apply_forces_to_cross_terms!, apply_forces!, apply_boundary_conditions!, linear_operator, update_diagonal_blocks!, update_offdiagonal_blocks!

function partial_update_state_dependents(sim_dd, s, m::SimulationModel, dt, g_forces; time = time)
    # Do nothing
end

function partial_update_linearized_system!(sim_dd, s, m::SimulationModel)
    # Do nothing
end


function partial_update_state_dependents(sim_dd, storage, model::MultiModel, dt, forces; time = time)
    p = sim_dd.partition
    targets, sources, all_models = submodel_keys_asymmetric(p, model)
    # Update all equations except for main / reservoir
    for key in targets
        Jutul.update_equations!(storage[key], model[key], dt)
    end
    # Then update the cross terms
    Jutul.update_cross_terms!(storage, model, dt, targets = all_models, sources = all_models)
    # Apply cross terms + bc to other equations than those of the the main model
    Jutul.apply_forces!(storage, model, dt, forces, time = time, targets = targets)
    Jutul.apply_forces_to_cross_terms!(storage, model, dt, forces; time = time, targets = all_models, sources = all_models)
    Jutul.apply_boundary_conditions!(storage, model, targets = targets)
    # Jutul.apply_cross_terms!(storage, model, dt, targets = targets, sources = sources)
end


function partial_update_linearized_system!(sim_dd, storage, model::MultiModel)
    p = sim_dd.partition
    targets, sources, all_models = submodel_keys_asymmetric(p, model)
    update_diagonal_blocks!(storage, model, targets)
    update_offdiagonal_blocks!(storage, model, all_models, sources)
end

function remove_local_update_from_aspen_residual!(sim_dd, storage)
    # Subtract right hand side A*Δv from aspen residual to skip post-process.
    sys = storage.LinearizedSystem
    r = vec(sys.r_buffer)
    op = linear_operator(sys, skip_red = true)
    dvar = vec(sys.dx_buffer)
    update_dx_with_full_update!(sim_dd, storage)
    mul!(r, op, dvar, -1, true)
end

function update_dx_with_full_update!(sim_dd, storage)
    lsys_g = storage.LinearizedSystem
    positions = sim_dd.storage.primary_variable_positions
    buffers = sim_dd.storage.primary_variable_buffers
    n = length(positions)
    for i = 1:n
        update_dx_with_full_update_subdomain!(lsys_g, positions[i], buffers[i])
    end
end

function update_dx_with_full_update_subdomain!(lsys_g, positions, buffers)
    for i in eachindex(positions)
        pos = positions[i]
        buf = buffers[i]
        L = lsys_g[i, i]
        dx = L.dx_buffer
        zippy!(dx, pos, buf)
    end
end

function zippy!(r, pos, buf)
    for (p, v) in zip(pos, buf)
        @inbounds r[p] = v
    end
end

function submodel_keys_asymmetric(p, model)
    main_symbol = p.main_symbol
    all_models = Jutul.submodels_symbols(model)
    targets = filter(!isequal(main_symbol), all_models)
    sources = all_models
    return (targets, sources, all_models)
end
