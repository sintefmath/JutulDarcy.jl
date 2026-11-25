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
    # Apply cross terms + bc to other equations than those of the main model
    Jutul.apply_forces!(storage, model, dt, forces, time = time, targets = targets)
    Jutul.apply_forces_to_cross_terms!(storage, model, dt, forces; time = time, targets = all_models, sources = all_models)
    Jutul.apply_boundary_conditions!(storage, model, targets = targets)
    # Jutul.apply_cross_terms!(storage, model, dt, targets = targets, sources = sources)
end


function partial_update_linearized_system!(sim_dd, storage, model::MultiModel)
    p = sim_dd.partition
    targets, sources, all_models = submodel_keys_asymmetric(p, model)
    # Cross terms are added to the Jacobian (not assigned), so the entries must
    # be zeroed first. The reservoir block is filled in afterwards by the ASPEN
    # assembly.
    Jutul.zero_jacobian_entries!(storage.LinearizedSystem)
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
    # With MPI: get the local updates of the ghost cells from their owners
    exchange_ghost_increments!(sim_dd)
    mul!(r, op, dvar, -1, true)
end

"""
    exchange_ghost_increments!(sim_dd)

When running with MPI, fill the increments (`dx_buffer`) of the ghost cells with
the local updates computed by the processes that own them. Each process then
has the full increment needed for the coupling term of the ASPEN rows. Does
nothing for a serial simulation.
"""
function exchange_ghost_increments!(sim_dd)
    executor = sim_dd.simulator.executor
    if executor isa Jutul.PArrayExecutor
        exchange = get(executor.data, :distributed_increment_exchange, nothing)
        if isnothing(exchange)
            error("ASPEN with MPI needs the increment exchange set up by setup_reservoir_simulator_parray.")
        end
        exchange()
    end
    return sim_dd
end

"""
    finalize_distributed_aspen_system!(outer_sim)

When running with MPI, replace the rows of ghost cells by unit rows, as is done
after a regular Newton assembly. Does nothing for a serial simulation.
"""
function finalize_distributed_aspen_system!(outer_sim)
    executor = outer_sim.executor
    if executor isa Jutul.PArrayExecutor
        s = outer_sim.storage
        Jutul.post_update_linearized_system!(s.LinearizedSystem, executor, s, outer_sim.model)
    end
    return outer_sim
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
