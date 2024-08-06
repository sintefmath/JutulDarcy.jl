function global_aspen!(
        simulator,
        dt,
        all_forces,
        forces, config,
        iteration,
        solve_status,
        subreports;
        relaxation = 1.0,
        report = Jutul.setup_ministep_report(),
        sim_kwarg...
    )
    outer_sim = simulator.simulator
    storage = simulator.storage
    s = outer_sim.storage
    m = outer_sim.model
    rec = Jutul.progress_recorder(simulator)
    time = rec.recorder.time + dt
    t_aspen_jac = t_aspen_gres = t_aspen_lres = t_aspen_couple = t_storage = t_conv = 0.0
    e = errors = nothing
    converged = false
    done_enough_its = iteration > config[:min_nonlinear_iterations]
    check_fi = false
    lsys = s.LinearizedSystem
    check = config[:debug_checks]
    sub_sims = simulator.subdomain_simulators
    n = length(sub_sims)
    function subdomain_solved_in_single_timestep(i)
        return solve_status[i] == local_solved_in_single_step
    end
    if false
        @info "FI fallback."
        e, converged, report = Jutul.perform_step!(s, m, dt, forces, config; iteration = iteration, update_secondary = true, sim_kwarg...)
        converged = converged && done_enough_its
    else
        buffers = storage.primary_variable_buffers
        pos = storage.primary_variable_positions
        p = simulator.partition
        t_storage = @elapsed store_subdomain_primary_variables!(simulator, config, do_delta = true)
        solve_for_full_increment = config[:aspen_full_increment]
        check_fi = solve_for_full_increment
        if check_fi
            t_secondary, t_equations = update_state_dependents!(s, m, dt, forces, time = time, update_secondary = true)
            t_linsys = @elapsed Jutul.update_linearized_system!(s, m)
            t_conv = @elapsed converged, e, errors = check_global_convergence(s, m, config, iteration, dt)
            Jutul.extra_debug_output!(report, s, m, config, iteration, dt)
        else
            # Update models to have the right equations and cross terms (leaving out the reservoir model
            # since that will be done in a special manner due to ASPEN.
            t_secondary = 0.0
            t_equations = @elapsed partial_update_state_dependents(simulator, s, m, dt, forces, time = time)
            # Update some parts of the linear system (everything except the big reservoir-reservoir block)
            t_linsys = @elapsed partial_update_linearized_system!(simulator, s, m)
        end
        converged = converged && done_enough_its

        if !converged
            # Update the residuals D(p^k - p^{k+1}) in each subdomain
            @tic "aspen local res" t_aspen_lres = @elapsed for i in 1:n
                sub_sim = sub_sims[i]
                if subdomain_solved_in_single_timestep(i)
                    update_aspen_residual_local!(sub_sim, buffers[i], i, p)
                else
                    reassemble_newton_locally!(s, outer_sim, sub_sim, all_forces, dt, i, buffers[i], pos[i])
                end
            end
            # Compute fluxes etc at the boundary of the local subdomains
            @tic "aspen cross-flux" t_aspen_couple = @elapsed for i in 1:n
                update_aspen_coupling!(storage, s, sub_sims[i].model, sub_sims[i].storage.state, i)
                check && check_aspen_coupling!(sub_sims[i].storage, sub_sims[i].model, storage.boundary_discretizations[i], i)
            end
            # Update global Jacobians with all local fluxes and what exists in each subdomain
            @tic "linear system" t_aspen_jac = @elapsed begin
                update_aspen_jacobian_global!(storage, outer_sim, sub_sims, p)
            end
            # Now fill in the global right hand side
            # Set rhs to D(p^k - p^{k+1/2}) -> rhs is now solving for p^k to p^{k+1} update
            @tic "aspen global res" t_aspen_gres = @elapsed for i in 1:n
                sub_sim = sub_sims[i]
                if subdomain_solved_in_single_timestep(i)
                    update_aspen_residual_global!(s, sub_sim.storage.LinearizedSystem, pos[i], i, p)
                end
            end
            # Subtract A*(p^k - p^{k+1/2}) from r to solve for only p^{k+1/2} to p^k
            @tic "aspen global preprocess" if !solve_for_full_increment
                # TODO: Check that this doesn't mess up stuff when FI fallback is used
                remove_local_update_from_aspen_residual!(simulator, s)
            end
            for i in 1:n
                sub_sim = sub_sims[i]
                if !subdomain_solved_in_single_timestep(i)
                    set_global_residual_to_newton_residual_for_subdomain!(s.LinearizedSystem, sub_sim, storage, i)
                end
            end
            t_conv = @elapsed if !check_fi
                @tic "convergence" converged, e, errors = check_global_convergence(s, m, config, iteration, dt)
                Jutul.extra_debug_output!(report, s, m, config, iteration, dt)
            end
            converged = converged && done_enough_its
        end
        # Store timings etc
        report[:secondary_time] = t_secondary
        report[:equations_time] = t_equations + t_storage + t_aspen_couple
        report[:linear_system_time] = t_linsys + t_aspen_jac + t_aspen_gres + t_aspen_lres
        report[:converged] = converged
        report[:errors] = errors
        report[:convergence_time] = t_conv
        if !converged
            lsolve = config[:linear_solver]
            check = config[:safe_mode]
            update_for_global_preconditioner!(simulator, lsolve, dt)
            lsys = s.LinearizedSystem
            t_solve = @elapsed begin
                @tic "linear solve" (lsolve_ok, n_iter, rep_lsolve) = Jutul.linear_solve!(lsys, lsolve, m, s, dt, rec)
            end
            if solve_for_full_increment
                t_aspen_lres = @elapsed @batch minbatch=10 for i in 1:n
                    if subdomain_solved_in_single_timestep(i)
                        subtract_primary_deltas!(simulator.simulator, buffers[i], pos[i])
                    end
                end
            end
            t_update = @elapsed @tic "primary variables" report[:update] = Jutul.update_primary_variables!(s, m, check = check; relaxation = relaxation)

            report[:linear_solver] = rep_lsolve
            report[:linear_iterations] = n_iter
            report[:linear_solve_time] = t_solve
            report[:update_time] = t_update
            report[:solved] = true
        else
            report[:solved] = false
        end
    end

    return e, converged, report
end

function update_for_global_preconditioner!(simulator, lsolve, dt)
    if isa(lsolve, GenericKrylov) && isa(lsolve.preconditioner, CPRPreconditioner) && lsolve.preconditioner.strategy == :true_impes
        # If we are using CPR with true impes we need the conservation law to be updated.
        rs = simulator.simulator.storage[:Reservoir]
        tm = rs.state.TotalMasses
        tm0 = rs.state0.TotalMasses
        acc = rs.equations.mass_conservation.accumulation.entries
        @. acc = (tm - tm0)/dt
    end
end

function check_global_convergence(s, m, config, iteration, dt)
    converged, e, errors = check_convergence(s, m, config, iteration = iteration, dt = dt, extra_out = true)
    # Act like the FI solver and produce output as requested
    il = config[:info_level]
    if il > 1
        Jutul.get_convergence_table(errors, il, iteration, config)
    end
    if iteration <= config[:min_nonlinear_iterations]
        converged = false
    end
    return (converged, e, errors)
end

function reassemble_newton_locally!(s, outer_sim, sub_sim, forces, dt, i, buffers = nothing, pos = nothing)
    # Update all variables from global since the boundary could have changed.
    # Secondary variables only depend on primary variables and are always in
    # sync in the local solves. If a time step fails, they are reset to the
    # values at the start.
    update_subdomain_from_global(outer_sim, sub_sim, i, current = true, transfer_secondary = true)
    s_i = sub_sim.storage
    m_i = sub_sim.model
    f_i = local_forces(forces, i)

    update_subdomain_from_global(outer_sim, sub_sim, i, current = false, transfer_secondary = true)
    Jutul.update_state_dependents!(s_i, m_i, dt, f_i, update_secondary = false)
    Jutul.update_linearized_system!(s_i, m_i)
    if !isnothing(buffers)
        for buffer in buffers
            @. buffer = 0.0
        end
        for (i, p) in enumerate(pos)
            dx_buf = s.LinearizedSystem[i, i].dx_buffer
            @inbounds for k in p
                dx_buf[k] = 0
            end
        end
    end
end

function update_aspen_residual_global!(storage_g, lsys, positions, p_i, p)
    lsys_g = storage_g.LinearizedSystem
    for i in eachindex(positions)
        L = lsys[i, i]
        pos = positions[i]
        buf = L.r_buffer
        r = lsys_g[i, i].r_buffer
        if Jutul.is_cell_major(L.matrix_layout)
            bz = Jutul.block_size(L)
            pos = equation_major_to_block_major_view(pos, bz)
        end
        zippy_set!(r, pos, buf)
    end
end

function zippy_set!(r, pos, buf)
    for (p, v) in zip(pos, buf)
        @inbounds r[p] = v
    end
end

function zippy_plus!(r, pos, buf)
    for (p, v) in zip(pos, buf)
        @inbounds r[p] += v
    end
end

function zippy_minus!(r, pos, buf)
    for (p, v) in zip(pos, buf)
        @inbounds r[p] -= v
    end
end

function zippy_swap_minus!(r, pos, buf)
    for (p, v) in zip(pos, buf)
        @inbounds r[p] = v - r[p]
    end
end

function subtract_primary_deltas!(sim, buffers, positions)
    storage_g = sim.storage
    lsys_g = storage_g.LinearizedSystem
    for i in eachindex(positions)
        L = lsys_g[i, i]
        pos = positions[i]
        buf = buffers[i]
        dx = lsys_g[i, i].dx_buffer

        if Jutul.is_cell_major(L.matrix_layout)
            # @info "Hey"
            # bz = Jutul.block_size(L)
            # pos = equation_major_to_block_major_view(pos, bz)
        end
        zippy_plus!(dx, pos, buf)
    end
end

function update_aspen_residual_local!(sim, buffer, i, partition)
    model = sim.model
    storage = sim.storage
    lsys = storage.LinearizedSystem
    update_aspen_residual_local!(lsys, buffer)
end

function update_aspen_residual_local!(lsys::Jutul.MultiLinearizedSystem, buffers)
    subsys = lsys.subsystems
    n = size(subsys, 1)
    # Do diagonals first (which zero out the buffers too)
    for i = 1:n
        update_aspen_residual_local!(subsys[i, i], buffers[i])
    end
    for j = 1:n
        sys = subsys[j, j]
        bz_d = element_block_size(sys.r)
        col_is_block = bz_d > 1
        for i = 1:n
            if i == j
                # Already treated and in-place.
                continue
            end
            r = subsys[i, i].r
            sys_ij = subsys[i, j]
            J = sys_ij
            op = linear_operator(sys_ij)
            bz = Jutul.block_size(sys_ij)
            buf = buffers[j]
            r_buf = reinterpret(Float64, r)
            if col_is_block && bz == 1
                # Column is in block_major but rhs is always row_major
                buf_b = equation_major_to_block_major_view(buf, bz_d)
            else
                # Straightforward, but recast for r
                buf_b = buf
            end
            # TODO: This collect could be preallocated. It is here due to
            # generic matmul fallback issues with reinterpreted arrays
            tmp = collect(r_buf)
            mul!(tmp, op, collect(buf_b), 1.0, 1.0)
            @. r_buf = tmp
        end
    end
    # @info "Inner: $t_in outer: $t_out"
end


element_block_size(r::AbstractArray{T}) where T<:Float64 = 1
element_block_size(r::AbstractArray{T}) where T = length(T)


function set_global_residual_to_newton_residual_for_subdomain!(lsys_g, sub_sim, storage, i)
    s_i = sub_sim.storage
    lsys_l = s_i.LinearizedSystem
    positions = storage.primary_variable_positions[i]
    for j in eachindex(positions)
        pos = positions[j]
        l = lsys_l[j, j]
        r = l.r_buffer
        if Jutul.is_cell_major(l.matrix_layout)
            bz = Jutul.block_size(l)
            pos = equation_major_to_block_major_view(pos, bz)
        end
        L = lsys_g[j, j]
        R = L.r_buffer
        zippy_set!(R, pos, r)
    end
end

function update_aspen_residual_local!(lsys::LinearizedSystem, buffers::Vector{Vector{Float64}})
    update_aspen_residual_local!(lsys, only(buffers))
end

function update_aspen_residual_local!(lsys::LinearizedSystem, buffer)
    J = lsys.jac
    r = lsys.r
    dx = lsys.dx_buffer
    update_aspen_residual_local_inner!(r, J, buffer, dx)
end

function update_aspen_residual_local_inner!(r0, J, buffer, dx_buf)
    # r = D * (u_{k+1/2} - u_k)
    # println("Inner $(eltype(r))")
    b, r, is_copy = mul_buffer(r0, buffer, dx_buf)
    mul!(r, J, b, 1, false)
    if is_copy
        @. r0 = r
    end
    return nothing
end

function mul_buffer(r::AbstractVector{T}, buffer, dx_buf) where T<:AbstractFloat
    return (buffer, r, false)
end

function mul_buffer(r::AbstractVector{T}, buffer, dx_buf) where T
    # Missing block-ordering
    bz = length(T)
    b_v = equation_major_to_block_major_view(buffer, bz)
    for i in eachindex(b_v)
        @inbounds dx_buf[i] = b_v[i]
    end
    b = reinterpret(T, vec(dx_buf))
    return (b::AbstractVector{T}, collect(r), true)
end

function postprocess_aspen_dx!(storage, s, m, submodel, substorage, positions, buffer, p_i, p)
    lsys = s.LinearizedSystem
    ng = length(buffer)
    for g = 1:ng
        L = lsys[g, g]
        dx = L.dx_buffer
        pos = positions[g]::Vector{Int64}
        buf = buffer[g]::Vector{Float64}

        zipcrement!(dx, pos, buf)
    end
end

function zipcrement!(dx, pos, buf)
    for (p, v) in zip(pos, buf)
        @inbounds dx[p] += v
    end
end

function equation_major_to_block_major_view(a, block_size)
    @views x = reshape(reshape(vec(a), :, block_size)', :)
    return x
end

function block_major_to_equation_major_view(a, block_size)
    @views x = reshape(reshape(vec(a), block_size, :)', :)
    return x
end

# Primary variable buffers
primary_variable_buffer(storage, model) = get_pvar_buffer(storage.LinearizedSystem)

function get_pvar_buffer(lsys::LinearizedSystem, is_wrapped = false)
    n = length(lsys.r_buffer)
    if is_wrapped
        o = zeros(n)
    else
        o = [zeros(n)]
    end
    return o
end

function get_pvar_buffer(lsys::Jutul.MultiLinearizedSystem)
    sub = lsys.subsystems
    n = size(sub, 1)
    return map(i -> get_pvar_buffer(sub[i, i], true), 1:n)
end

# Positions for primary variable buffers
function primary_variable_positions(outer_sim, sim, buffers, i, partition)
    primary_variable_positions(outer_sim.model, sim.model, sim.storage, buffers, i, partition)
end

function primary_variable_positions(model_g::MultiModel, model::MultiModel, storage, buffers::Vector{Vector{Float64}}, i, partition)
    sym = Jutul.submodels_symbols(model_g)
    groups = model_g.groups
    single_group = isnothing(groups)
    if single_group
        groups = ones(Int64, length(sym))
    end
    models_g = model_g.models
    models = model.models
    # Some models may just have the first group (e.g. skipping wells)
    ng = min(maximum(groups), length(buffers))
    pos_g = Vector{Vector{Int64}}()
    for g in 1:ng
        # For each group and corresponding linear system:
        # - loop over the models in global model

        # buffer = buffers[g]
        # offset = 0
        subs = findall(isequal(g), groups)
        ctx = models_g[sym[first(subs)]].context

        # Determine layout of group by first entry
        layout = Jutul.matrix_layout(ctx)
        cell_major = Jutul.is_cell_major(layout)
        nb = length(buffers[g])
        # Should put this somewhere.
        pos = zeros(Int64, nb)

        global_offset = 0
        # local_offset = 0
        current_pos = 1
        for ix in subs
            s = sym[ix]
            m_g = models_g[s]
            if haskey(models, s)
                m = models[s]
            else
                m = nothing
            end
            global_offset, current_pos = store_model_primary_variables!(pos, m_g, m, cell_major, global_offset, current_pos)
        end
        # @info "Group $g:" cell_major nb buffers pos
        push!(pos_g, pos)
        @assert nb + 1 == current_pos "$nb+1 should be equal to $current_pos"
    end
    return pos_g
end

function primary_variable_positions(model_g, model, storage, buffers, i, partition)
    buf = buffers[1]
    nb = length(buf)
    pos = zeros(Int64, nb)
    layout = Jutul.matrix_layout(model.context)
    cell_major = Jutul.is_cell_major(layout)
    store_model_primary_variables!(pos, model_g, model, cell_major, 0, 1)
    return [pos]
end

function store_model_primary_variables!(pos, m_g, m, cell_major, global_offset, current_pos)
    if isnothing(m)
        d = nothing
        primary = nothing
    else
        d = m.domain
        primary = m.primary_variables
    end
    d_g = m_g.domain
    primary_g = m_g.primary_variables

    if cell_major
        # Assumes one unit with primary variables in the block backend
        for u in Jutul.get_primary_variable_ordered_entities(m)
            n_global = Jutul.count_entities(d_g, u)
            # Total number of partials in entity
            block_size = Jutul.number_of_partials_per_entity(m, u)
            # Local -> global map
            active = Jutul.active_entities(d, u)
            inner_offset = 0
            for (pkey, p) in primary_g
                if !haskey(primary, pkey)
                    continue
                end
                if u != Jutul.associated_entity(p)
                    continue
                end
                shift = global_offset + inner_offset
                np = Jutul.degrees_of_freedom_per_entity(m, p)
                for j = 1:np
                    for (local_i, c_i) in enumerate(active)
                        global_i = Jutul.global_cell(c_i, global_map(d))
                        # DD primary storage always follows equation_major style
                        next = (global_i-1)*block_size + j + shift
                        pos[current_pos] = next
                        current_pos += 1
                    end
                end
                inner_offset += np
            end
            global_offset += n_global*block_size
        end
    else
        for (pkey, p) in primary_g
            u = Jutul.associated_entity(p)
            # How much space does this variable take up in the global ordering?
            # We need this even if the variable is not present in the submodel we 
            # are looking at to get correct indices.
            n_global = Jutul.count_active_entities(d_g, u)
            np = Jutul.degrees_of_freedom_per_entity(m_g, p)
            width_global = n_global*np
            if !isnothing(primary) && haskey(primary, pkey)
                n_local = Jutul.count_active_entities(d, u)
                # Local -> global map
                active = Jutul.active_entities(d, u)
                for j = 1:np
                    for c_i in active
                        global_i = Jutul.global_cell(c_i, global_map(d))
                        pos[current_pos] = (j-1)*n_global + global_offset + global_i
                        current_pos += 1
                    end
                end
            else
                n_local = 0
            end
            global_offset += width_global
        end
    end
    return (global_offset, current_pos)
end


# Storage of primary variables in ASPEN
function store_subdomain_primary_variables!(simulator, config; kwarg...)
    s = simulator.storage
    buffers = s.primary_variable_buffers
    if haskey(s, :black_oil_primary_buffers)
        extra = s[:black_oil_primary_buffers]
    else
        extra = nothing
    end
    for i in eachindex(buffers)
        if isnothing(extra)
            e = nothing
        else
            e = extra[i]
        end
        store_primary_variables!(simulator.subdomain_simulators[i], buffers[i], e; kwarg...)
    end
end

function store_primary_variables!(sim::JutulSimulator, buffer, extra; kwarg...)
    store_primary_variables!(sim.model, sim.storage, buffer, extra; kwarg...)
end

# store_primary_variables!(model::SimulationModel, sim, buffer::Vector{Vector{Float64}}; kwarg...) = store_primary_variables!(model, sim, buffer[1]; kwarg...)

function store_primary_variables!(model::SimulationModel, storage, buffer::Vector{Vector{Float64}}, extra_buffer; kwarg...)
    @assert length(buffer) == 1
    # @info "Writing to buffer 1"
    if isnothing(extra_buffer)
        e = nothing
    else
        e = extra_buffer[1]
    end
    store_primary_variables!(model::SimulationModel, storage, buffer[1], e; kwarg...)
end

function store_primary_variables!(model::SimulationModel, storage, buffer::Vector{Float64}, extra; offset = 0, do_delta = false)
    pvars = storage.primary_variables
    pvar_base = model.primary_variables
    state = storage.state
    # @info "starting buffer write" size(buffer)
    offset0 = offset
    # @info keys(pvar_base) typeof(pvar_base)
    for k in keys(pvar_base)
        # @info "Variable $k"
        pb = pvar_base[k]
        pvar = pvars[k]
        base_scaling = @something Jutul.variable_scale(pb) 1.0
        scale = 1.0/base_scaling
        # TODO: Fix allocations here.
        a = active_for_transfer(model.domain, pvar, pb)
        np = transfer_pvar!(buffer, model, state, a, do_delta, pvar, pb, extra, scale, offset)
        offset += np
    end
    return offset - offset0
end

active_for_transfer(d::Jutul.DiscretizedDomain, pvar, pb) = Jutul.active_entities(d, Jutul.associated_entity(pb))
active_for_transfer(d, pvar, pb) = 1:number_of_entities(pvar)

function transfer_pvar!(buffer, model, state, a, do_delta, pvar::AbstractVector, pb, extra, scale, offset)
    n = length(a)
    # @info "$(typeof(pb))" n
    @inbounds for i = 1:n
        v = value(pvar[a[i]])*scale
        pos = i + offset
        # @info "Writing to $pos <vec> (Δ = $do_delta)"
        if do_delta
            @inbounds v = buffer[pos] - v
        end
        buffer[pos] = v
    end
    return n
end

function transfer_pvar!(buffer, model, state, a, do_delta, pvar::AbstractMatrix, pb, extra, scale, offset)
    n = length(a)
    m = Jutul.degrees_of_freedom_per_entity(model, pb)::Int
    @inbounds for j = 1:m
        for i = 1:n
            # Assume that it is the m first entries that are primary
            v = value(pvar[j, a[i]])*scale
            pos = (j-1)*n + i + offset
            # @info "Writing to $pos <mat> (Δ = $do_delta)"
            if do_delta
                @inbounds v = buffer[pos] - v
            end
            buffer[pos] = v
        end
    end
    return n*m
end

function transfer_pvar!(buffer, model, state, a, do_delta, pvar::AbstractVector, pb::JutulDarcy.BlackOilUnknown, extra, scale, offset)
    @assert !isnothing(extra)
    @assert scale == 1.0
    n = length(a)
    sys = model.system
    vapoil = JutulDarcy.has_vapoil(sys)
    disgas = JutulDarcy.has_disgas(sys)

    check = false
    @inline function check_val(old, new, i, sym)
        if check
            new_val = new.val
            if abs(old - new_val) > 1e-8
                @error "Bad value in $i for $sym:" old new new_val
                error("Stopping.")
            end
        end
    end
    # TODO: Not 100% robust
    if haskey(model.secondary_variables, :Rs)
        boreg = model[:Rs].regions
    elseif haskey(model.secondary_variables, :Rv)
        boreg = model[:Rv].regions
    else
        boreg = nothing
    end

    @inbounds for (i, cell) in enumerate(a)
        p = value(state.Pressure[cell])
        if haskey(state, :ImmiscibleSaturation)
            swi = value(state.ImmiscibleSaturation[cell])
        else
            swi = 0.0
        end
        reg_i = JutulDarcy.region(boreg, cell)
        v = pvar[cell]
        phases = v.phases_present
        x = value(v.val)
        if disgas
            if phases == JutulDarcy.OilOnly
                rs = x
            elseif phases == JutulDarcy.OilAndGas
                rstab = JutulDarcy.table_by_region(model.system.rs_max, reg_i)
                rs = rstab(p)
            else
                # Oil only
                rs = 0.0
            end
        else
            rs = 0.0
        end
        if vapoil
            if phases == JutulDarcy.GasOnly
                rv = x
            elseif phases == JutulDarcy.OilAndGas
                rvtab = table_by_region(model.system.rs_max, reg_i)
                rv = rvtab(p)
            else
                # Gas only
                rv = 0.0
            end
        else
            rv = 0.0
        end
        if phases == JutulDarcy.OilAndGas
            # Note: sg as a primary variable is really
            # Sg/(1 - sw), i.e. a value for which 1.0 
            # means that there is no oil.
            sg = x
        elseif phases == JutulDarcy.GasOnly
            sg = 1.0 - swi
        else
            sg = 0.0
        end
        pos = i + offset
        if do_delta
            # Note: Reverse order due to Jacobian etc
            # @inbounds v = buffer[pos] - v
            if phases == JutulDarcy.OilOnly && disgas
                rs0 = extra[1, i]
                check_val(rs, v, i, :rs)
                Δ = rs0 - rs
            elseif phases == JutulDarcy.GasOnly && vapoil
                rv0 = extra[2, i]
                Δ = rv0 - rv
            else # Two-phase
                sg0 = extra[3, i]
                check_val(sg, v, i, :sg)
                Δ = sg0 - sg
            end
            buffer[pos] = Δ
        else
            phases = v.phases_present
            if phases == JutulDarcy.OilOnly && disgas
                check_val(rs, v, i, :rs)
            elseif phases == JutulDarcy.GasOnly && vapoil
                check_val(rv, v, i, :rv)
            else # Two-phase
                check_val(sg, v, i, :sg)
            end
            # rs, rv, sg
            extra[1, i] = rs
            extra[2, i] = rv
            extra[3, i] = sg
            buffer[pos] = 0
        end
    end
    return n
end

function store_primary_variables!(model::MultiModel, storage, buffers, extra; offset = 0, do_delta = false)
    @assert offset == 0
    sym = Jutul.submodel_symbols(model)
    groups = model.groups
    if isnothing(groups)
        groups = ones(Int64, length(sym))
        # buffers = [buffers]
    end
    ng = maximum(groups)
    for g in 1:ng
        # @info "Writing to buffer $g"
        buffer = buffers[g]
        offset = 0
        subs = findall(isequal(g), groups)
        for (i, ix) in enumerate(subs)
            s = sym[ix]
            if isnothing(extra)
                e = nothing
            else
                e = extra[s]
            end
            # @info "$s"
            n = store_primary_variables!(model.models[s], storage[s], buffer, e, offset = offset, do_delta = do_delta)
            offset += n
        end
    end
    return offset
end

function write_system(lsys, name)
    # myprint(io, i, j, v) = write(io, @sprintf("%d, %d: %1.4g\n", ctr, v))
    @info "Writing"
    if isa(lsys, Jutul.MultiLinearizedSystem)
        for i = 1:2
            for j = 1:2
                open("$(name)_system_$(i)_$(j)", "w") do io
                    L = lsys.subsystems[i, j]
                    jac = L.jac
                    I, J, V = findnz(jac)
                    for (row, col, v) in zip(I, J, V)
                        write(io, "$row, $col: ")
                        write(io, "$v\n")
                    end
                end
            end
            open("$(name)_rhs_$(i)", "w") do io
                L = lsys.subsystems[i, i]
                r = L.r
                for (ix, v) in enumerate(r)
                    write(io, "$ix: $v\n")
                end
            end
        end
    else
        error()
    end
    error("Wrote system")
end

function print_sys_order_indep(lsys, suffix)
    myprint(io, ctr, v) = write(io, @sprintf("%d: %1.4g\n", ctr, v))
    if isa(lsys, Jutul.MultiLinearizedSystem)
        open("block$suffix.txt", "w") do io
            for i = 1:2
                L = lsys.subsystems[i, i]
                if use_r
                    r = L.r
                    rb = L.r_buffer
                else
                    r = L.dx
                    rb = L.dx_buffer
                end
                n = length(r)
                if i == 1
                    for j = 1:2
                        for c = 1:n
                            v = rb[(c-1)*2 + j]
                            myprint(io, ctr, v)
                            ctr += 1
                        end
                    end
                else
                    # @info i > 2
                    for v in r
                        myprint(io, ctr, v)
                        ctr += 1
                    end
                end
            end
        end
    else
        open("sparse$suffix.txt", "w") do io
            if use_r
                r = lsys.r
            else
                r = lsys.dx
            end
            for v in r
                myprint(io, ctr, v)
                ctr += 1
            end
        end
    end
end


function print_rhs_order_indep(lsys, use_r, suffix)
    if use_r
        @info "Writing r"
    else
        @info "Writing dx"
    end
    ctr = 1
    # myprint(io, ctr, v) = write(io, "$ctr: $v\n")
    myprint(io, ctr, v) = write(io, @sprintf("%d: %1.4g\n", ctr, v))

    # myprint(io, ctr, v) = println("$ctr: $v")

    if isa(lsys, Jutul.MultiLinearizedSystem)
        open("block$suffix.txt", "w") do io
            for i = 1:2
                L = lsys.subsystems[i, i]
                if use_r
                    r = L.r
                    rb = L.r_buffer
                else
                    r = L.dx
                    rb = L.dx_buffer
                end
                n = length(r)
                if i == 1
                    for j = 1:2
                        for c = 1:n
                            v = rb[(c-1)*2 + j]
                            myprint(io, ctr, v)
                            ctr += 1
                        end
                    end
                else
                    # @info i > 2
                    for v in r
                        myprint(io, ctr, v)
                        ctr += 1
                    end
                end
            end
        end
    else
        open("sparse$suffix.txt", "w") do io
            if use_r
                r = lsys.r
            else
                r = lsys.dx
            end
            for v in r
                myprint(io, ctr, v)
                ctr += 1
            end
        end
    end
end

function check_aspen_assembly(storage, sub_sims, p)
    sys = storage.LinearizedSystem[1, 1]
    part = p.partition[:Reservoir].partition
    for (i, s) in enumerate(sub_sims)
        ix = findall(part .== i)
        g_part = sys.jac[ix, ix]
        l_part = s.storage.LinearizedSystem[1, 1].jac
        d = g_part - l_part
        # @info "Block $i" g_part l_part g_part - l_part
        @assert length(nonzeros(d)) == 0
    end
end
