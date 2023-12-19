include("processing.jl")
function ijk_to_linear(i, j, k, dims)
    nx, ny, nz = dims
    return (k-1)*nx*ny + (j-1)*nx + i
end

function ij_to_linear(i, j, dims)
    nx, ny = dims
    return (j-1)*nx + i
end

function linear_to_ijk(ix, dims)
    nx, ny, nz = dims
    linear_index = ix
    x = mod(linear_index - 1, nx) + 1
    y = mod((linear_index - x) ÷ nx, ny) + 1
    leftover = (linear_index - x - (y-1)*nx)
    z = (leftover ÷ (nx*ny)) + 1
    return (x, y, z)
end

function get_line(coord, i, j, nx, ny)
    ix = ijk_to_linear(i, j, 1, (nx, ny, 1))
    T = SVector{3, Float64}
    x1 = T(coord[ix, 1:3])
    x2 = T(coord[ix, 4:end])

    return (x1, x2)
end

# function interp_coord(line::Tuple, z)
#     p0, p1 = line
#     return interp_coord(p0, p1, z)
# end

function interp_coord(line::NamedTuple, z)
    if line.equal_points
        x = line.x1
        T = eltype(x)
        pt = SVector{3, T}(x[1], x[2], z)
    else
        pt = interp_coord(line.x1, line.x2, z)
    end
    return pt
end

function interp_coord(p0::SVector{3, T}, p1::SVector{3, T}, z::T) where T<:Real
    z0 = p0[3]
    z1 = p1[3]
    if z0 ≈ z1 
        # Coinciding corner points. Just return the point and hope the pillar is
        # inactive.
        interp_pt = p0
    else
        weight = (z - z0)/(z1 - z0)
        interp_pt = p0 .+ weight.*(p1 .- p0)
        @assert interp_pt[3] ≈ z "expected $z was $(interp_pt[3]) != $z"
    end
    return interp_pt
end

function corner_index(cell, corner::NTuple, dims)
    # Cell 1 [corner1, corner2], cell 2 [corner1, corner2], ..., cell n [corner1, corner2] repeated for nx in top plane
    # Cell 1 [corner3, corner4], cell 2 [corner3, corner4], ..., cell n [corner3, corner4]
    # Cell 1 [corner5, corner6], cell 2 [corner5, corner6], ..., cell n [corner5, corner6]
    # Cell 1 [corner7, corner8], cell 2 [corner7, corner8], ..., cell n [corner7, corner8]

    if cell isa Int
        i, j, k = linear_to_ijk(cell, dims)
    else
        i, j, k = cell
        cell = ijk_to_linear(i, j, k, cartdims)
    end
    nx, ny, nz = dims
    # offset = 2*prod(dims)
    # near/far, left/right, up/down
    if true
        i_is_upper, j_is_upper, k_is_upper = corner

        j_is_upper, i_is_upper, k_is_upper = corner
        cell_i_offset = 2*(i-1)
        cell_j_offset = 4*nx*(j-1)
        cell_k_offset = 8*nx*ny*(k-1)
        # @info "big offsets" cell_i_offset cell_j_offset cell_k_offset

        cell_offset = cell_i_offset + cell_j_offset + cell_k_offset

        i_offset = i_is_upper+1
        j_offset = j_is_upper*2*nx
        k_offset = k_is_upper*4*nx*ny
        # @info "small offsets" i_offset j_offset k_offset

        ijk = i_offset + j_offset + k_offset
        return cell_offset + ijk

        # x = i_offset + i_upper
        # y = j_offset + j_upper*2*nx
        # z = k_offset + k_upper*4*nx*ny
        # return x + y + z
        # return k*k_offset + j*j_offset + i_offset + 


        # i_pos = i_offset + x + 1
        # j_pos = 

        # offset = i_offset + j_offset + k_offset

        # return planar_offset + top_bottom_offset + y_offset + x_offset + minor
    
    end
    error()
    # if corner == (0, 0, 0)
    #     major = 0
    #     minor = 1
    # elseif corner == (1, 0, 0)
    #     major = 0
    #     minor = 2
    # elseif corner == (0, 1, 0)
    #     major = 1
    #     minor = 1
    # elseif corner == (1, 1, 0)
    #     major = 1
    #     minor = 2
    # elseif corner == (0, 0, 1)
    #     major = 2
    #     minor = 1
    # elseif corner == (0, 1, 1)
    #     major = 2
    #     minor = 2
    # elseif corner == (1, 0, 1)
    #     major = 3
    #     minor = 1
    # elseif corner == (1, 1, 1)
    #     major = 3
    #     minor = 2
    # else
    #     error("Unsupported $corner_type")
    # end
    # if corner == (0, 0, 0)
    #     major = 0
    #     minor = 1
    # elseif corner == (1, 0, 0)
    #     major = 1
    #     minor = 1
    # elseif corner == (0, 1, 0)
    #     major = 0
    #     minor = 2
    # elseif corner == (1, 1, 0)
    #     major = 1
    #     minor = 2
    # elseif corner == (0, 0, 1)
    #     major = 2
    #     minor = 1
    # elseif corner == (0, 1, 1)
    #     major = 2
    #     minor = 2
    # elseif corner == (1, 0, 1)
    #     major = 3
    #     minor = 1
    # elseif corner == (1, 1, 1)
    #     major = 3
    #     minor = 2
    # else
    #     error("Unsupported $corner_type")
    # end
    planar_offset = 8*nx*ny*(k-1)
    top_bottom_offset = major*2*nx*ny
    y_offset = 2*nx*(j-1)
    x_offset = 2*(i-1)
    return planar_offset + top_bottom_offset + y_offset + x_offset + minor
end


function cpgrid_primitives(coord, zcorn, cartdims; actnum = missing)
    # Add all lines that have at least one active neighbor
    coord = reshape(coord, 6, :)'
    nx, ny, nz = cartdims
    if ismissing(actnum)
        actnum = Array{Bool, 3}(undef, nx, ny, nz)
        @. actnum = true
    end
    # remapped_indices = findall(vec(actnum))
    nactive = sum(vec(actnum))
    remapped_indices = Vector{Int}(undef, nx*ny*nz)
    tmp = vec(actnum)
    active_cell_indices = findall(isequal(1), tmp)
    @. remapped_indices[tmp] = 1:nactive

    nlinex = nx+1
    nliney = ny+1
    @assert nliney*nlinex == size(coord, 1)

    function generate_line(p1, p2)
        return (
            z = Vector{Float64}(),
            cells = Vector{Int}(),
            cellpos = Vector{Int}(),
            nodes = Vector{Int}(),
            x1 = SVector{3, Float64}(p1),
            x2 = SVector{3, Float64}(p2),
            equal_points = p1 ≈ p2
            )
    end
    # active_lines = BitArray(undef, nlinex, nliney)
    x1, x2 = get_line(coord, 1, 1, nlinex, nliney)
    line0 = generate_line(x1, x2)
    function boundary_index(i, j, is_top)
        if is_top
            layer_offset = -(2*nx*ny*nz)
        else
            layer_offset = -(nx*ny*nz)
        end
        return layer_offset #- ij_to_linear(i, j, cartdims[1:2])
    end

    function cell_index(i, j, k)
        ix = ijk_to_linear(i, j, k, cartdims)
        if actnum[i, j, k]
            cell = remapped_indices[ix]
        else
            cell = -ix
            @assert cell <= 0
        end
        return cell
    end

    linear_line_ix(i, j) = ij_to_linear(i, j, (nlinex, nliney))
    lines = Matrix{typeof(line0)}(undef, nlinex, nliney)
    for i in 1:nlinex
        for j in 1:nliney
            p1, p2 = get_line(coord, i, j, nlinex, nliney)
            lines[i, j] = generate_line(p1, p2)
        end
    end
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                ix = ijk_to_linear(i, j, k, cartdims)
                active_cell_index = cell_index(i, j, k)
                for I1 in (0, 1)
                    for I2 in (0, 1)
                        L = lines[i + I2, j + I1]
                        for I3 in (0, 1)
                            zcorn_ix = corner_index(ix, (I1, I2, I3), cartdims)
                            c = zcorn[zcorn_ix]
                            # Note reversed indices, this is a bit of a mess
                            push!(L.z, c)
                            push!(L.cells, active_cell_index)
                        end
                    end
                end
            end
        end
    end
    # Add fake boundary nodes with corresponding cells
    for i in 1:(nx+1)
        for j in 1:(ny+1)
            L = lines[i, j]
            # Top layer
            if L.cells[1] > 0
                # t = boundary_index(i, j, true)
                # pushfirst!(L.z, L.z[1] - 1.0, L.z[1])
                # pushfirst!(L.cells, t, t)
            end
            # Bottom layer
            if L.cells[end] > 0
                # b = boundary_index(i, j, false)
                # push!(L.z, L.z[end], L.z[end] + 1.0)
                # push!(L.cells, b, b)
            end
        end
    end

    # Process lines and merge similar nodes
    nodes, lines_active = process_lines!(lines)

    # The four lines making up each column
    column_lines = Vector{NTuple{4, Int64}}()
    # for i = 1:nx
    #     for j = 1:ny
    #         ll = linear_line_ix(i, j)
    #         rl = linear_line_ix(i+1, j)
    #         lr = linear_line_ix(i, j+1)
    #         rr = linear_line_ix(i+1, j+1)

    #         ll_act = lines_active[ll]
    #         rl_act = lines_active[rl]
    #         lr_act = lines_active[lr]
    #         rr_act = lines_active[rr]

    #         if ll_act && rl_act && lr_act && rr_act
    #             push!(column_lines, (ll, rl, rr, lr))
    #         end
    #     end
    # end
    # Tag columns as active or inactive
    active_columns = Matrix{Bool}(undef, nx, ny)
    for i in 1:nx
        for j in 1:ny
            is_active = false
            for k in 1:nz
                is_active = is_active || actnum[i, j, k]
            end
            active_columns[i, j] = is_active
        end
    end
    # Generate the columns with cell lists
    make_column(i, j) = (cells = Int[], i = i, j = j)
    cT = typeof(make_column(1, 1))
    col_counter = 1
    columns = Vector{cT}()
    column_indices = zeros(Int, nx, ny)
    for i in 1:nx
        for j in 1:ny
            if active_columns[i, j]
                col = make_column(i, j)
                prev = boundary_index(i, j, true) 
                # push!(col.cells, prev)
                for k in 1:nz
                    cell = cell_index(i, j, k)
                    if cell != prev
                        push!(col.cells, cell)
                    end
                    prev = cell
                end
                # Put a boundary at the end
                # push!(col.cells, boundary_index(i, j, false))
                push!(columns, col)
                column_indices[i, j] = col_counter
                col_counter += 1

                ll = linear_line_ix(i, j)
                rl = linear_line_ix(i+1, j)
                lr = linear_line_ix(i, j+1)
                rr = linear_line_ix(i+1, j+1)
                push!(column_lines, (ll, rl, rr, lr))
            end
        end
    end
    ncol = length(columns)
    ncoll = length(column_lines)
    @assert ncol == ncoll "Mismatch in columns ($ncol) and column lines ($ncoll)"

    function get_edge(i, j, t)
        if t == :right
            p1 = linear_line_ix(i+1, j+1)
            p2 = linear_line_ix(i+1, j)
        elseif t == :left
            p1 = linear_line_ix(i, j)
            p2 = linear_line_ix(i, j+1)
        elseif t == :upper
            p1 = linear_line_ix(i, j+1)
            p2 = linear_line_ix(i+1, j+1)
        else
            @assert t == :lower
            p1 = linear_line_ix(i, j)
            p2 = linear_line_ix(i+1, j)
        end
        return (p1, p2)
    end

    function get_boundary_edge(self, i, j, t)
        return (column = self, pillars = get_edge(i, j, t))
    end

    function get_interior_edge(c1, c2, i, j, t)
        return (columns = (c1, c2), pillars = get_edge(i, j, t))
    end

    tmp = get_boundary_edge(1, 1, 1, :left)
    column_boundary = Vector{typeof(tmp)}()

    tmp = get_interior_edge(1, 1, 1, 1, :left)
    column_neighbors = Vector{typeof(tmp)}()

    for i in 1:nx
        for j in 1:ny
            if active_columns[i, j]
                self = column_indices[i, j]
                if i < nx && active_columns[i+1, j]
                    other = column_indices[i+1, j]
                    e = get_interior_edge(self, other, i, j, :right)
                    push!(column_neighbors, e)
                else
                    # Add right edge to boundary
                    e = get_boundary_edge(self, i, j, :right)
                    push!(column_boundary, e)
                end
                if i == 1
                    e = get_boundary_edge(column_indices[i, j], i, j, :left)
                    push!(column_boundary, e)
                end
                if j < ny && active_columns[i, j+1]
                    other = column_indices[i, j+1]
                    e = get_interior_edge(self, other, i, j, :upper)
                    push!(column_neighbors, e)
                else
                    e = get_boundary_edge(self, i, j, :upper)
                    push!(column_boundary, e)
                end
                if j == 1
                    e = get_boundary_edge(column_indices[i, j], i, j, :lower)
                    push!(column_boundary, e)
                end
            else
                if i < nx && active_columns[i+1, j]
                    e = get_boundary_edge(column_indices[i+1, j], i+1, j, :left)
                    push!(column_boundary, e)
                end
                if j < ny && active_columns[i, j+1]
                    e = get_boundary_edge(column_indices[i, j+1], i, j+1, :lower)
                    push!(column_boundary, e)
                end
            end
        end
    end

    return (
        lines = lines,
        lines_active = lines_active,
        column_neighbors = column_neighbors,
        column_boundary = column_boundary,
        column_lines = column_lines,
        columns = columns,
        active = active_cell_indices,
        nodes = nodes,
        cartdims = cartdims
    )
end

function process_lines!(lines)
    nodes = Vector{SVector{3, Float64}}()
    active_lines = BitVector(undef, length(lines))
    node_counter = 1
    for (line_ix, line) in enumerate(lines)
        z = line.z
        active = length(z) > 0 && !all(x -> x <= 0, line.cells)
        active_lines[line_ix] = active
        if active
            p = sortperm(z)
            @. line.z = z[p]
            @. line.cells = line.cells[p]
            pos = line.cellpos
            push!(pos, 1)

            counter = 1
            for i in 2:length(z)
                if z[i] ≈ z[i-1]
                    counter += 1
                else
                    push!(pos, pos[end] + counter)
                    counter = 1
                end
            end
            push!(pos, pos[end] + counter)
            # Sort each set of cells
            for i in 1:(length(pos)-1)
                ci = view(line.cells, pos[i]:(pos[i+1]-1))
                sort!(ci)
            end
            ix = pos[1:end-1]
            unique_z = line.z[ix]
            # Put back the unique points only
            resize!(line.z, 0)
            for z_i in unique_z
                push!(line.z, z_i)
                push!(line.nodes, node_counter)
                node_counter += 1
                new_node = interp_coord(line, z_i)
                push!(nodes, new_node)
            end
        end
    end

    return (nodes, active_lines)
end

function grid_from_primitives(primitives)
    (;
        lines,
        lines_active,
        # column_neighbors,
        column_lines,
        columns,
        active,
        nodes,
        cartdims
    ) = primitives
    # Faces mapping to nodes
    faces = Vector{Int}()
    face_pos = [1]
    faceno = 1

    cell_is_boundary(x) = x < 1
    # Boundary faces mapping to nodes
    boundary_faces = Vector{Int}()
    boundary_face_pos = [1]
    boundary_faceno = 1

    # Mapping from cell to faces
    cell_faces = Vector{Vector{Int}}()
    # Mapping from cell to boundary faces
    cell_boundary_faces = Vector{Vector{Int}}()

    for c in eachindex(active)
        push!(cell_faces, Vector{Int}())
        push!(cell_boundary_faces, Vector{Int}())
    end
    face_neighbors = Vector{Tuple{Int, Int}}()
    boundary_cells = Vector{Int}()

    nx, ny, nz = cartdims

    nlinex = nx+1
    nliney = ny+1

    extra_node_lookup = Dict()

    function add_face_from_nodes!(V, Vpos, nodes)
        n_global_nodes = length(primitives.nodes)
        n_local_nodes = length(nodes)
        @assert n_local_nodes > 2
        for n in nodes
            @assert n <= n_global_nodes
            @assert n > 0
            push!(V, n)
        end
        push!(Vpos, length(nodes) + Vpos[end])
    end

    function insert_boundary_face!(prev_cell, cell, nodes)
        orient = cell_is_boundary(prev_cell) && !cell_is_boundary(cell)
        @assert orient || (cell_is_boundary(cell) && !cell_is_boundary(prev_cell)) "cell pair $((cell, prev_cell)) is not on boundary"
        if orient
            self = cell
        else
            self = prev_cell
            nodes = reverse(nodes)
        end
        add_face_from_nodes!(boundary_faces, boundary_face_pos, nodes)
        push!(cell_boundary_faces[self], boundary_faceno)
        push!(boundary_cells, self)
        boundary_faceno += 1
    end

    function insert_interior_face!(prev_cell, cell, nodes)
        @assert cell > 0
        @assert prev_cell > 0
        @assert prev_cell != cell
        add_face_from_nodes!(faces, face_pos, nodes)
        # Note order here.
        push!(face_neighbors, (prev_cell, cell))
        push!(cell_faces[cell], faceno)
        push!(cell_faces[prev_cell], faceno)
        faceno += 1
    end

    function insert_face!(prev_cell, cell, nodes; is_boundary)
        if is_boundary
            insert_boundary_face!(prev_cell, cell, nodes)
        else
            insert_interior_face!(prev_cell, cell, nodes)
        end
    end

    # Horizontal faces (top/bottom and faces along column)
    node_buffer = Int[]
    sizehint!(node_buffer, 10)
    for (cl, col) in zip(column_lines, columns)
        number_of_cells_in_column = length(col.cells)
        current_column_lines = map(l -> lines[l], cl)
        for (i, cell) in enumerate(col.cells)
            if cell_is_boundary(cell)
                continue
            end
            if i == 1
                prev = 0
            else
                prev = col.cells[i-1]
            end
            if i == number_of_cells_in_column
                next = 0
            else
                next = col.cells[i+1]
            end
            top_is_boundary = cell_is_boundary(prev)
            bottom_is_boundary = cell_is_boundary(next)
            cell_bnds = map(l -> find_cell_bounds(cell, l), current_column_lines)
            for is_top in (true, false)
                # TODO: c1/c2 definition will have to be modified to pick the
                # right normal vector, check this later.
                if is_top
                    if !top_is_boundary
                        # Avoid adding interior faces twice.
                        continue
                    end
                    is_bnd = top_is_boundary
                    F = first
                    c1 = prev
                    c2 = cell
                else
                    is_bnd = bottom_is_boundary
                    F = last
                    c1 = cell
                    c2 = next
                end
                # Index into pillars
                node_in_pillar_indices = map(F, cell_bnds)
                # Then find the global node indices
                node_indices = map((line, i) -> line.nodes[i], current_column_lines, node_in_pillar_indices)
                insert_face!(c1, c2, node_indices, is_boundary = is_bnd)
            end
        end
    end
    # primitives.column_boundary
    for is_bnd in [true, false]
        if is_bnd
            col_neighbors = primitives.column_boundary
        else
            col_neighbors = primitives.column_neighbors
        end
        for (cols, pillars) in col_neighbors
            # Get the pair of lines we are processing
            p1, p2 = pillars
            l1 = lines[p1]
            l2 = lines[p2]
            if length(cols) == 1
                a = b = only(cols)
            else
                a, b = cols
            end

            col_a = columns[a]
            col_b = columns[b]

            cell_pairs, overlaps = JutulDarcy.traverse_column_pair(col_a, col_b, l1, l2)
            int_pairs, int_overlaps, bnd_pairs, bnd_overlaps = split_overlaps_into_interior_and_boundary(cell_pairs, overlaps)

            F_interior = (l, r, node_indices) -> insert_face!(l, r, node_indices, is_boundary = false)
            F_bnd = (l, r, node_indices) -> insert_face!(l, r, node_indices, is_boundary = true)

            if is_bnd
                # We are dealing with a boundary column, everything is boundary
                add_vertical_cells_from_overlaps!(extra_node_lookup, F_bnd, nodes, int_pairs, int_overlaps, l1, l2)
            else
                add_vertical_cells_from_overlaps!(extra_node_lookup, F_interior, nodes, int_pairs, int_overlaps, l1, l2)
                add_vertical_cells_from_overlaps!(extra_node_lookup, F_bnd, nodes, bnd_pairs, bnd_overlaps, l1, l2)
            end
        end
    end

    function convert_to_flat(v)
        flat_vals = Int[]
        flat_pos = Int[1]
        for cf in v
            for face in cf
                push!(flat_vals, face)
            end
            push!(flat_pos, flat_pos[end]+length(cf))
        end
        return (flat_vals, flat_pos)
    end

    c2f, c2f_pos = convert_to_flat(cell_faces)
    b2f, b2f_pos = convert_to_flat(cell_boundary_faces)

    return UnstructuredMesh(
        c2f,
        c2f_pos,
        b2f,
        b2f_pos,
        faces,
        face_pos,
        boundary_faces,
        boundary_face_pos,
        primitives.nodes,
        face_neighbors,
        boundary_cells;
        structure = CartesianIndex(cartdims[1], cartdims[2], cartdims[3]),
        cell_map = primitives.active
    )
end

function next_cell_line_ptr!(line_ptr, self_lines, next_cell)
    for (i, line_i) in enumerate(self_lines)
        cells_i = line_i.cells
        ptr_i = line_ptr[i]
        cellpos = line_i.cellpos
        current_i = cellpos[ptr_i]
        nc = length(cells_i)
        while current_i <= nc && cells_i[current_i] != next_cell
            current_i += 1
        end
        next_ptr = -1
        for j in ptr_i:(length(cellpos)-1)
            if current_i >= cellpos[j] && current_i < cellpos[j+1]
                next_ptr = j
                break
            end
        end
        if next_ptr == -1
            @info "?!" line_i cells_i next_cell
            error("This should not occur.")
        end
        pos = cellpos[next_ptr]:(cellpos[next_ptr+1]-1)
        @assert next_cell in cells_i[pos]
        line_ptr[i] = next_ptr
    end
end