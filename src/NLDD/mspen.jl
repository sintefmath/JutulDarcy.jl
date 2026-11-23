# MSPEN: Multiplicative version of ASPEN.
#
# The local stage is a Gauss-Seidel sweep where subdomain i is solved with the
# subdomains earlier in the sweep already updated. Write the block Jacobian
# with rows of subdomain i evaluated where subdomain i was solved as
# J = D + 𝓛 + 𝓤, where 𝓛 is the coupling to subdomains earlier in the sweep and
# 𝓤 the coupling to later subdomains. The Jacobian of the preconditioned system
# is then (D + 𝓛)⁻¹J (compared to D⁻¹J for ASPEN) and the global system has the
# same matrix as ASPEN. The right hand side differs from the ASPEN one by 𝓛b,
# where b is the local update, which is added here. The result -𝓤b is also the
# linearized residual at the state after the sweep.

function add_mspen_lower_coupling!(simulator, storage_g, order; check = false)
    sub_sims = simulator.subdomain_simulators
    n = length(sub_sims)
    sweep_rank = fill(n + 1, n)
    for (k, i) in enumerate(order)
        sweep_rank[i] = k
    end
    owner = mspen_cell_owner(simulator)
    lsys = storage_g.LinearizedSystem[1, 1]
    r = vec(lsys.r_buffer)
    dx = vec(lsys.dx_buffer)
    nc = length(owner)
    bz = length(r) ÷ nc
    @assert bz*nc == length(r)
    cell_major = Jutul.is_cell_major(lsys.matrix_layout)
    if check
        r_before = copy(r)
    end
    model_g = mspen_reservoir_model(simulator.simulator.model)
    laws = aspen_conservation_laws(model_g)
    for i in 1:n
        bnd = simulator.storage.boundary_discretizations[i]
        M = global_map(mspen_reservoir_model(sub_sims[i].model).domain)
        conn_data = bnd.mass_flow.conn_data
        for (key, bflux) in pairs(bnd.flux)
            add_mspen_lower_coupling_law!(r, dx, Jutul.get_entries(bflux), conn_data, M, owner, sweep_rank, i, laws[key], nc, bz, cell_major)
        end
    end
    if check
        check_mspen_lower_coupling(simulator, lsys, owner, sweep_rank, r .- r_before, dx, nc, bz, cell_major)
    end
    return r
end

function add_mspen_lower_coupling_law!(r, dx, ∂v, conn_data, M, owner, sweep_rank, i, equation_offset, nc, bz, cell_major)
    rank_i = sweep_rank[i]
    ne = size(∂v, 1)
    for (j, cd) in enumerate(conn_data)
        source = Jutul.global_cell(cd.self, M)
        if sweep_rank[owner[source]] < rank_i
            target = Jutul.global_cell(cd.other, M)
            for e in 1:ne
                # The Jacobian entries are the negative partials of the
                # boundary flux (see update_from_boundary_fluxes!)
                ∂q = ∂v[e, j].partials
                v = 0.0
                for d in eachindex(∂q)
                    v -= ∂q[d]*dx[mspen_linear_index(source, d, nc, bz, cell_major)]
                end
                r[mspen_row_index(target, e, equation_offset, nc, bz, cell_major)] += v
            end
        end
    end
end

function mspen_linear_index(cell, k, nc, bz, cell_major)
    if cell_major
        ix = (cell - 1)*bz + k
    else
        ix = (k - 1)*nc + cell
    end
    return ix
end

function mspen_row_index(cell, e, equation_offset, nc, bz, cell_major)
    # The equation offset is the number of rows of the preceding equations
    # (counted over all cells), as used by find_jac_position.
    if cell_major
        ix = (cell - 1)*bz + equation_offset ÷ nc + e
    else
        ix = equation_offset + (e - 1)*nc + cell
    end
    return ix
end

mspen_reservoir_model(m::MultiModel) = m.models.Reservoir
mspen_reservoir_model(m) = m

function mspen_cell_owner(simulator)
    s = simulator.storage
    if !haskey(s, :mspen_cell_owner)
        sub_sims = simulator.subdomain_simulators
        nc = number_of_cells(mspen_reservoir_model(simulator.simulator.model).domain)
        owner = zeros(Int, nc)
        for (i, sim) in enumerate(sub_sims)
            D = mspen_reservoir_model(sim.model).domain
            M = global_map(D)
            for c in Jutul.active_entities(D, M, Cells())
                owner[Jutul.global_cell(c, M)] = i
            end
        end
        @assert all(x -> x > 0, owner) "All cells must belong to a subdomain."
        s[:mspen_cell_owner] = owner
    end
    return s[:mspen_cell_owner]
end

function check_mspen_lower_coupling(simulator, lsys, owner, sweep_rank, added, dx, nc, bz, cell_major)
    # Compute 𝓛b by masked products with the assembled global Jacobian
    op = linear_operator(lsys)
    n = length(simulator.subdomain_simulators)
    expected = zeros(length(added))
    cell_of_index(ix) = cell_major ? div(ix - 1, bz) + 1 : mod(ix - 1, nc) + 1
    for j in 1:n
        b_j = zeros(length(dx))
        for ix in eachindex(b_j)
            if owner[cell_of_index(ix)] == j
                b_j[ix] = dx[ix]
            end
        end
        y = op*b_j
        for ix in eachindex(y)
            i = owner[cell_of_index(ix)]
            if sweep_rank[j] < sweep_rank[i]
                expected[ix] += y[ix]
            end
        end
    end
    err = norm(added - expected, Inf)/max(norm(expected, Inf), 1e-20)
    jutul_message("MSPEN", "Relative error in lower coupling: $err (|𝓛b| = $(norm(expected, Inf)))")
    return err
end
