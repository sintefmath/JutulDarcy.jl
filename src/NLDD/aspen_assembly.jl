function update_aspen_jacobian_global!(storage, outer_sim, sub_sims, mp)
    # We only need to deal with the reservoir due to assumptions on partitioning for other things
    s = outer_sim.storage
    m = outer_sim.model

    n = length(sub_sims)
    @batch minbatch = 10 for i = 1:n
        sim = sub_sims[i]
        bflux = storage.boundary_discretizations[i]
        update_aspen_jacobian_global!(bflux, m, s, sim.model, sim.storage, i)
    end
end

function update_aspen_jacobian_global!(bflux, model_g, storage_g, model, storage, i)
    # We only need to deal with the reservoir due to assumptions on partitioning for other things
    lsys_g = storage_g.LinearizedSystem[1, 1]
    lsys = storage.LinearizedSystem[1, 1]
    get_resmodel(m::MultiModel, s) = (m.models.Reservoir, s.Reservoir)
    get_resmodel(m, s) = (m, s)

    rm_g, rs_g = get_resmodel(model_g, storage_g)
    rm, rs = get_resmodel(model, storage)

    update_aspen_jacobian_global_reservoir!(lsys_g, lsys, bflux, rm_g, rs_g, rm, rs, i)
end

function update_aspen_jacobian_global_reservoir!(lsys_g, lsys, bflux, model_g, storage_g, model, storage, block_no)
    # If the global Jacobian has already been assembled, we can skip most of this. Thus, if the outer code
    # is ever made smarter, the following must be set to false.
    nz_g = lsys_g.jac_buffer
    nz = lsys.jac_buffer
    # Global version
    law_g = model_g.equations[:mass_conservation]
    law_storage_g = storage_g.equations.mass_conservation
    # Local version
    law = model.equations[:mass_conservation]
    law_storage = storage.equations.mass_conservation
    acc_g = law_storage_g.accumulation
    acc = law_storage.accumulation

    _, ne, np = Jutul.ad_dims(acc_g)
    D = model.domain
    M = global_map(D)
    active_cells = Jutul.active_entities(D, M, Cells())

    fill_diagonals!(nz_g, nz, active_cells, M, acc_g, acc, np, ne)

    # Fill in the partials of fluxes in the interior
    flux_g = law_storage_g.half_face_flux_cells
    fpos_g = flux_g.jacobian_positions

    flux = law_storage.half_face_flux_cells
    fentries = flux.entries
    fpos = flux.jacobian_positions

    conn_pos = law.flow_discretization.conn_pos
    conn_data = law.flow_discretization.conn_data

    conn_pos_g = law_g.flow_discretization.conn_pos
    conn_data_g = law_g.flow_discretization.conn_data

    fill_interior_fluxes(nz_g, flux, fentries, active_cells, fpos, flux_g, fpos_g, M, conn_pos, conn_data, conn_pos_g, conn_data_g, ne, np)
    # Finally, add the lagged fluxes on the boundary of the interior
    bdisc = bflux.mass_flow
    bflux_cache = bflux.flux
    ∂v = Jutul.get_entries(bflux_cache)
    ∂pos = bflux_cache.jacobian_positions
    nb = size(∂v, 2)
    update_from_boundary_fluxes!(nz_g, bdisc, bflux_cache, ∂pos, ∂v, nb, ne, np)
end

function fill_diagonals!(nz_g, nz, active_cells, M, acc_g, acc, np, ne)
    @inbounds for (i, c_i) in enumerate(active_cells)
        g_i = Jutul.global_cell(c_i, M)
        for d in 1:np
            @inbounds for e in 1:ne
                # look up global pos
                g_pos = Jutul.get_jacobian_pos(acc_g, g_i, e, d)
                # look up local pos
                pos = Jutul.get_jacobian_pos(acc, i, e, d)
                # insert into global buffer from local
                @inbounds V = nz[pos]
                @inbounds nz_g[g_pos] = V
            end
        end
    end
end

function update_from_boundary_fluxes!(nz_g, bdisc, bflux_cache, ∂pos, ∂v, nb, ne, np)
    for j in 1:nb
        for e in 1:ne
            @inbounds v = -∂v[e, j]
            @inbounds for d in 1:np
                face_pos_global = Jutul.get_jacobian_pos(bflux_cache, j, e, d, ∂pos)
                nz_g[face_pos_global] = v.partials[d]
            end
        end
    end
end

function fill_interior_fluxes(nz_g, flux, fentries, active_cells, fpos, flux_g, fpos_g, M, conn_pos, conn_data, conn_pos_g, conn_data_g, ne, np)
    @inbounds for (i, c_i) in enumerate(active_cells)
        g_i = Jutul.global_cell(c_i, M)
        start = conn_pos[i]
        stop = conn_pos[i+1]
        n_f = stop - start
        @inbounds for j in 0:(n_f-1)
            f = conn_data[start+j].face
            f_g = M.faces[f]
            # @info "Local cell $i ($c_i) for $block_no accessing $(i+j) to get local face $f global $f_g"
            # Find the same global face in the global ordering of half faces for the corresponding global cell
            hf_pos_g = 0
            @inbounds for k = conn_pos_g[g_i]:conn_pos_g[g_i+1]-1
                cd_g = conn_data_g[k]
                if cd_g.face == f_g
                    hf_pos_g = k
                end
            end
            if hf_pos_g == 0
                # Something is wrong.
                error("Did not find face $f_g for global cell $g_i in subdomain $block_no local cell $c_i")
            end
            # Local position is much easier
            hf_pos = start + j
            done_with_face = false
            # Loop over all inner faces in the subdomain and insert into the global Jacobian.
            # We do not have to do double work filling in the diagonal since that was already
            # copied over above.
            for e in 1:ne
                @inbounds q = -Jutul.get_entry(flux, hf_pos, e)
                @inbounds for d in 1:np
                    face_pos_local = Jutul.get_jacobian_pos(flux, hf_pos, e, d, fpos)
                    if face_pos_local == 0
                        # Exterior boundary - skip along
                        done_with_face = true
                        break
                    end
                    face_pos_global = Jutul.get_jacobian_pos(flux_g, hf_pos_g, e, d, fpos_g)
                    @inbounds ∂q = q.partials[d]
                    @inbounds nz_g[face_pos_global] = ∂q
                end
                if done_with_face
                    break
                end
            end
        end
    end
end