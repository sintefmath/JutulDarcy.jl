"""
    reservoir_transmissibility(d::DataDomain)
    reservoir_transmissibility(d::DataDomain; version = :xyz)

Special transmissibility function for reservoir simulation that handles
additional complexity present in industry grade models such as fault
multipliers, net-to-gross etc. The input should be a `DataDomain` instance
returned from [`reservoir_domain`](@ref)

The keyword argument `version` can be `:xyz` for permeability tensors that are
aligned with coordinate directions or `:ijk` to interpreted the permeability as
a diagonal tensor aligned with the logical grid directions. The latter choice is
only meaningful for a diagonal tensor.
"""
function reservoir_transmissibility(d::DataDomain;
        version = :xyz,
        projection = missing
    )
    nf = number_of_faces(d)
    has_nnc = haskey(d, :nnc)
    if has_nnc
        nnc = d[:nnc]
        nnc::NonNeighboringConnections
        num_nnc = length(nnc.trans_flow)
    else
        nnc = missing
        num_nnc = 0
    end

    function apply_ntg!(T_hf, ntg, facepos, face_is_vertical)
        for (c, ntg) in enumerate(ntg)
            if ntg isa AbstractFloat && ntg ≈ 1.0
                continue
            end
            for fp in facepos[c]:(facepos[c+1]-1)
                if face_is_vertical[faces[fp]]
                    T_hf[fp] *= ntg
                end
            end
        end
    end
    function fig_negative_trans!(T_hf, faceno)
        neg_count = 0
        for (i, T_hf_i) in enumerate(T_hf)
            if faceno[i] > nf - num_nnc
                continue
            end
            neg_count += T_hf_i < 0
            T_hf[i] = abs(T_hf_i)
        end
        # We only warn for significant amounts of negative transmissibilities, since
        # a few negative values is normal for reservoir grids.
        if neg_count > 0.1*length(T_hf)
            tran_tot = length(T_hf)
            perc = round(100*neg_count/tran_tot, digits = 2)
            jutul_message("Transmissibility", "Replaced $neg_count negative half-transmissibilities (out of $tran_tot, $perc%) with their absolute value.")
        end
        return T_hf
    end
    function fix_bad_trans!(T_hf, faceno)
        bad_count = 0
        for (i, T_hf_i) in enumerate(T_hf)
            if faceno[i] > nf - num_nnc
                continue
            end
            if T_hf_i isa AbstractFloat && !isfinite(T_hf_i)
                bad_count += 1
                T_hf[i] = 0.0
            end
        end
        if bad_count > 0
            tran_tot = length(T_hf)
            perc = round(100*bad_count/tran_tot, digits = 2)
            jutul_message("Transmissibility", "Replaced $bad_count non-finite half-transmissibilities (out of $tran_tot, $perc%) with zero.")
        end
        return T_hf
    end
    g = physical_representation(d)
    N = d[:neighbors]
    nc = number_of_cells(d)
    faces, facepos = get_facepos(N, nc)
    facesigns = Jutul.get_facesigns(N, faces, facepos, nc)

    if version == :ijk
        face_dir = get_ijk_face_dir(g, N)
    else
        face_dir = missing
    end
    if ismissing(projection)
        projection = :normal
        project_face_centroids = version == :ijk
    end
    if project_face_centroids
        half_face_centroids = project_half_face_centroids(d, projection)
    else
        half_face_centroids = missing
    end
    T_hf = compute_half_face_trans(
        d[:cell_centroids],
        d[:face_centroids],
        d[:normals],
        d[:areas],
        d[:permeability],
        faces, facepos, facesigns,
        version = version,
        face_dir = face_dir,
        half_face_centroids = half_face_centroids
    )
    nf = number_of_faces(d)
    fig_negative_trans!(T_hf, faces)
    fix_bad_trans!(T_hf, faces)
    if haskey(d, :net_to_gross)
        # Net to gross applies to vertical trans only
        otag = get_mesh_entity_tag(g, Faces(), :orientation, throw = false)
        if !ismissing(otag)
            # Use tags if provided
            face_is_vertical = map(1:nf) do face
                return mesh_entity_has_tag(g, Faces(), :orientation, :vertical, face)
            end
        elseif g isa CartesianMesh
            # Cartesian mesh is simple
            k_index = map(c -> cell_ijk(g, c), 1:nc)
            face_is_vertical = map(1:nf) do face
                l, r = N[:, face]
                return k_index[l] == k_index[r]
            end
        else
            # Fall back to normals
            normals = d[:normals]
            face_is_vertical = map(1:nf) do face
                nx, ny, nz = normals[:, face]
                return abs(nz) < max(abs(nx), abs(ny))
            end
        end
        apply_ntg!(T_hf, d[:net_to_gross], facepos, face_is_vertical)
    end
    T = compute_face_trans(T_hf, N)
    if haskey(d, :transmissibility_multiplier, Faces())
        tm = d[:transmissibility_multiplier]
        @. T *= tm
    end
    if haskey(d, :transmissibility_override, Faces())
        for (f, v) in enumerate(d[:transmissibility_override])
            if isfinite(v)
                T[f] = v
            end
        end
    end
    if haskey(d, :numerical_aquifers)
        aquifers = d[:numerical_aquifers]
        bnd_areas = d[:boundary_areas]
        bnd_centroids = d[:boundary_centroids]
        cell_centroids = d[:cell_centroids]
        D = size(cell_centroids, 1)
        point_t = SVector{D, eltype(cell_centroids)}

        if haskey(d, :net_to_gross)
            ntg = d[:net_to_gross]
        else
            ntg = ones(nc)
        end
        T, num_aquifer_faces = set_aquifer_transmissibilities!(
            T,
            g, d[:permeability], ntg,
            aquifers,
            reinterpret(point_t, cell_centroids),
            reinterpret(point_t, bnd_centroids),
            g.boundary_faces.neighbors,
            bnd_areas
        )
    else
        num_aquifer_faces = 0
    end

    if has_nnc
        # NNC come at the end.
        offset = nf - num_nnc
        for (i, T_nnc) in enumerate(nnc.trans_flow)
            T[i + offset] = T_nnc
        end
    end
    return T
end

function set_aquifer_transmissibilities!(T, mesh, perm, ntg, aquifers, cell_centroids, bnd_centroids, bnd_neighbors, bnd_areas)
    num_aquifer_faces = 0
    # Connections to the reservoir
    for (aq_id, aquifer) in pairs(aquifers)
        aqprm = aquifer.aquifer_cells[1]
        aquifer_cell = aqprm.cell
        R = aqprm.length/2.0
        for (bface, face, opt, tmult) in zip(
                aquifer.boundary_faces,
                aquifer.added_faces,
                aquifer.trans_option,
                aquifer.boundary_transmult
            )
            area_reservoir = bnd_areas[bface]
            reservoir_cell = bnd_neighbors[bface]
            dist = norm(bnd_centroids[bface] - cell_centroids[reservoir_cell])
            num_aquifer_faces += 1
            is_vertical = mesh_entity_has_tag(mesh, BoundaryFaces(), :orientation, :vertical, bface)

            if mesh_entity_has_tag(mesh, BoundaryFaces(), :ijk_orientation, :j, bface)
                dir = 2
            elseif mesh_entity_has_tag(mesh, BoundaryFaces(), :ijk_orientation, :k, bface)
                dir = 3
            else
                dir = 1
            end
            if is_vertical
                ntg_face = ntg[reservoir_cell]
            else
                ntg_face = 1.0
            end
            T_reservoir = perm[dir, reservoir_cell]*area_reservoir*ntg_face/dist

            if opt == 0
                area_aquifer = aqprm.area
            else
                @assert opt == 1 "Option for aquifer transmissibility expected to be 1 or 0, was $opt"
                area_aquifer = area_reservoir
            end
            T_aquifer = area_aquifer*aqprm.permeability/R
            effective_trans = tmult/(1.0/T_reservoir + 1.0/T_aquifer)
            if isfinite(effective_trans)
                T[face] = effective_trans
            else
                @error "Non-finite aquifer transmissibility for numerical aquifer $aq_id, setting to zero" T_aquifer T_reservoir aqprm
                T[face] = 0.0
            end
        end
    end
    # Aquifer internal connections
    for (aq_id, aquifer) in pairs(aquifers)
        aqprms = aquifer.aquifer_cells
        @assert length(aquifer.aquifer_faces) == length(aqprms)-1
        for (i, face) in enumerate(aquifer.aquifer_faces)
            num_aquifer_faces += 1
            curr = aqprms[i]
            next = aqprms[i+1]
            T_c = curr.area*curr.permeability/(curr.length/2.0)
            T_n = next.area*next.permeability/(next.length/2.0)
            T[face] = 1.0/(1.0/T_c + 1.0/T_n)
        end
    end
    return (T, num_aquifer_faces)
end

function get_ijk_face_dir(g, N)
    nf = number_of_faces(g)
    face_dir = ones(Int, nf)
    ijk_c = map(i -> cell_ijk(g, i), 1:number_of_cells(g))
    gdim = dim(g)
    for i in 1:nf
        l, r = N[:, i]
        ijk_l = ijk_c[l]
        ijk_r = ijk_c[r]
        for j in 1:gdim
            if ijk_l[j] != ijk_r[j]
                face_dir[i] = j
                break
            end
        end
    end
    return face_dir
end

function project_half_face_centroids(
        face_centroids::AbstractVector,
        cell_centroids::AbstractVector,
        normals::AbstractVector,
        neighbors::AbstractMatrix{Int},
        hfc::AbstractVector{Int},
        hff::AbstractVector{Int},
        projection::Symbol = :normal
    )
    projection == :normal || throw(ArgumentError("Only :normal projection is supported"))
    nhf = length(hfc)
    out = sizehint!(similar(face_centroids, 0), nhf)
    for (cell, face) in zip(hfc, hff)
        hf = face_centroids[face]
        hc = cell_centroids[cell]
        sgn = neighbors[1, face] == cell ? 1.0 : -1.0
        normal = sgn*normals[face]

        proj_point = hf + dot(hc - hf, normal) * normal
        push!(out, proj_point)
    end
    return out
end

function project_half_face_centroids(d::DataDomain, projection = :normal)
    face_centroids = vec(reinterpret(SVector{3, Float64}, d[:face_centroids]))
    cell_centroids = vec(reinterpret(SVector{3, Float64}, d[:cell_centroids]))
    normals = vec(reinterpret(SVector{3, Float64}, d[:normals]))
    neighbors = d[:neighbors]

    hfc = d[:half_face_cells]
    hff = d[:half_face_faces]
    return project_half_face_centroids(face_centroids, cell_centroids, normals, neighbors, hfc, hff, projection)
end
