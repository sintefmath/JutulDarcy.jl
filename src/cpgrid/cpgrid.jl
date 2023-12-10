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
    if corner == (0, 0, 0)
        major = 0
        minor = 1
    elseif corner == (0, 0, 1)
        major = 2
        minor = 1
    elseif corner == (0, 1, 1)
        major = 2
        minor = 2
    elseif corner == (0, 1, 0)
        major = 0
        minor = 2
    elseif corner == (1, 0, 0)
        major = 1
        minor = 1
    elseif corner == (1, 0, 1)
        major = 3
        minor = 1
    elseif corner == (1, 1, 1)
        major = 3
        minor = 2
    elseif corner == (1, 1, 0)
        major = 1
        minor = 2
    else
        error("Unsupported $corner_type")
    end
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
                get_corner(pos) = zcorn[corner_index(ix, pos, cartdims)]
                for I1 in (0, 1)
                    for I2 in (0, 1)
                        L = lines[i + I2, j + I1]
                        for I3 in (0, 1)
                            c = get_corner((I1, I2, I3))
                            # Note reversed indices, this is a bit of a mess
                            push!(L.z, c)
                            push!(L.cells, cell_index(i, j, k))
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
                t = boundary_index(i, j, true)
                pushfirst!(L.z, L.z[1] - 1.0, L.z[1])
                pushfirst!(L.cells, t, t)
            end
            # Bottom layer
            if L.cells[end] > 0
                b = boundary_index(i, j, false)
                push!(L.z, L.z[end], L.z[end] + 1.0)
                push!(L.cells, b, b)
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
    make_column(i, j) = (cells = Int[typemin(Int)], i = i, j = j)
    cT = typeof(make_column(1, 1))
    col_counter = 1
    columns = Vector{cT}()
    column_indices = zeros(Int, nx, ny)
    for i in 1:nx
        for j in 1:ny
            if active_columns[i, j]
                col = make_column(i, j)
                prev = boundary_index(i, j, true) 
                push!(col.cells, prev)
                for k in 1:nz
                    cell = cell_index(i, j, k)
                    if cell != prev
                        push!(col.cells, cell)
                    end
                    prev = cell
                end
                # Put a boundary at the end
                push!(col.cells, boundary_index(i, j, false))
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
        column_neighbors,
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

    function add_face_from_nodes!(V, Vpos, nodes)
        for n in nodes
            push!(V, n)
        end
        push!(Vpos, length(nodes) + Vpos[end])
    end

    function insert_boundary_face!(prev_cell, cell, nodes)
        orient = cell_is_boundary(prev_cell) && !cell_is_boundary(cell)
        @assert orient || (cell_is_boundary(cell) && !cell_is_boundary(prev_cell))
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
    for (cl, col) in zip(column_lines, columns)
        # Traverse each column one by one and then figure out what cells are
        # connected, generate faces and push to the respective arrays
        ncorners = length(cl)
        line_ptr = ones(Int, ncorners)
        self_lines = map(i -> lines[i], cl)


        cells = col.cells
        @assert cell_is_boundary(cells[1])
        @assert cell_is_boundary(cells[end])
        cell = cells[2]
        # Deal with the top first
        nodes_first = map(x -> first(x.nodes), self_lines)
        n = length(cells)
        for i in 2:(n-1)
            prev_cell = cells[i-1]
            current_cell = cells[i]
            next_cell = cells[i+1]
            if cell_is_boundary(current_cell)
                # Keep going until we find actual cells
                continue
            end
            next_cell_line_ptr!(line_ptr, self_lines, current_cell)
            bottom_nodes = map((l, ix) -> l.nodes[ix+1], self_lines, line_ptr)
            if cell_is_boundary(prev_cell)
                # There was a gap or we were at the top, add boundary
                top_nodes = map((l, ix) -> l.nodes[ix], self_lines, line_ptr)
                insert_face!(prev_cell, current_cell, top_nodes; is_boundary = true)
            end
            # Next cell down in column might be a gap or at the bottom
            is_bnd = cell_is_boundary(next_cell)
            insert_face!(current_cell, next_cell, bottom_nodes; is_boundary = is_bnd)
        end
        # error()
    end
    # Vertical faces (between active columns)
    function initialize_cpair(pillars, cols, i)
        c = columns[cols[i]]
        @assert c.cells[1] <= 0
        return (lines[pillars[i]], c, c.cells[2])
    end
    
    function node_positions_cell_pair(cA, cB, line)
        # Find start and stop
        n = length(line.cellpos)-1
        get_pos(i) = line.cellpos[i]:(line.cellpos[i+1]-1)
        get_cells(i) = @view line.cells[get_pos(i)]
        start = missing
        for i in 1:n
            cells = get_cells(i)
            @info "$i: $cells"
            if cA in cells && cB in cells
                start = i
                break
            end
        end
        @assert !ismissing(start) "Did not find $((cA, cB)) in $(line.cells)"
        # @info "Found starting point" start
        stop = missing
        for i in n:-1:start
            cells = get_cells(i)
            # @info "??" cA cB cells i
            if cA in cells && cB in cells
                stop = i
                break
            end
        end
        @assert !ismissing(stop)
        # @info "Found ending point" stop n
        # pos = findfirst(isequal(cA), col.cells)
    
        if stop == n
            A_next = B_next = false
        else
            next = get_cells(stop+1)
            A_next = cA in next
            B_next = cB in next
        end
        @assert !(A_next && B_next) "Found $cA and $cB in $next at $(stop+1)?"
        # A in next, B in next
        # start and stop
        return (start, stop, A_next, B_next)
    end

    node_buf = zeros(Int, 10)
    for (cols, pillars) in column_neighbors
        lineA, colA, current_cellA = initialize_cpair(pillars, cols, 1)
        lineB, colB, current_cellB = initialize_cpair(pillars, cols, 2)
        ptrA = ptrB = 1
        colposA = colposB = 2
        # @info "?!" colA lineA

        nA = length(lineA.z)
        nB = length(lineB.z)
        it = 1
        # @warn "Starting" current_cellA current_cellB nA nB
        # @info "Current pillars" cols pillars
        # for i in 1:nA
        #     p = lineA.cellpos
        #     @error "Node $i:" lineA.cells[p[i]:(p[i+1]-1)]
        # end
        while true

            # Both functions should get both cell pairs
            startA, stopA, lineA_foundA, lineA_foundB = node_positions_cell_pair(current_cellA, current_cellB, lineA)
            startB, stopB, lineB_foundA, lineB_foundB = node_positions_cell_pair(current_cellA, current_cellB, lineB)

            # @info "Found?" foundA1 foundB1 ptrA ptrB next_ptrA next_ptrB current_cellA current_cellB

            # Insert new face
            bndA = cell_is_boundary(current_cellA)
            bndB = cell_is_boundary(current_cellB)
            if !bndA || !bndB
                Arange = startA:stopA
                Brange = startB:stopB
                resize!(node_buf, 0)
                for i in Arange
                    push!(node_buf, lineA.nodes[i])
                end
                for i in Iterators.reverse(Brange)
                    push!(node_buf, lineB.nodes[i])
                end
                # @info "Inserting face" primitives.nodes[node_buf] Arange Brange
                is_bnd = bndA || bndB
                insert_face!(current_cellA, current_cellB, node_buf; is_boundary = is_bnd)
            end
            # @info "??" lineA_foundA lineA_foundB lineB_foundA lineB_foundB
            # @warn "Iteration $it:" zA zB
            foundA1 = lineA_foundA || lineB_foundA
            if !foundA1
                if current_cellA == colA.cells[end]
                    break
                end
                @info "??" colA.cells current_cellA
                ok = false
                for i in eachindex(colA.cells)
                    if colA.cells[i] == current_cellA
                        current_cellA = colA.cells[i+1]
                        ok = true
                        break
                    end
                end
                if !ok
                    @info "Went bad" current_cellA colA
                    error("!!")
                end
            end
            foundB1 = lineA_foundB || lineB_foundB
            if !foundB1
                if current_cellB == colB.cells[end]
                    break
                end
                for i in eachindex(colB.cells)
                    if colB.cells[i] == current_cellB
                        current_cellB = colB.cells[i+1]
                        break
                    end
                end
            end
            @info "$it" current_cellA colA current_cellB colB foundB1 foundA1
            it += 1
        end

        # error()
    end
    # Treat boundary
    for (col, pillars) in primitives.column_boundary
        # lineA, lineB = pillars
        lineA, colA, current_cell = initialize_cpair(pillars, (col, col), 1)
        lineB, colB, current_cellB = initialize_cpair(pillars, (col, col), 2)
        @assert current_cell == current_cellB
        ptrA = ptrB = 1
        colpos = 2
        # @info "?!" colA lineA

        nA = length(lineA.z)
        nB = length(lineB.z)
        it = 1
        # @warn "Starting" current_cellA current_cellB nA nB
        # for i in 1:nA
        #     p = lineA.cellpos
        #     @error "Node $i:" lineA.cells[p[i]:(p[i+1]-1)]
        # end
        while ptrA < nA && ptrB < nB
            # next_ptrA, foundA1, foundB1 = seek_line(lineA, ptrA, current_cell, current_cell)
            # next_ptrB, foundA2, foundB2 = seek_line(lineB, ptrB, current_cell, current_cell)
            startA, stopA, foundA1, foundB1 = node_positions_cell_pair(current_cell, current_cell, lineA)
            startB, stopB, foundA2, foundB2 = node_positions_cell_pair(current_cell, current_cell, lineB)

            @assert !(foundA1 && foundA2 && foundB1 && foundB2)
            @assert foundA1 == foundA2
            @assert foundB1 == foundB2

            if !cell_is_boundary(current_cell)
                Arange = startA:stopA
                Brange = startB:stopB
                resize!(node_buf, 0)
                for i in Arange
                    push!(node_buf, lineA.nodes[i])
                end
                for i in Iterators.reverse(Brange)
                    push!(node_buf, lineB.nodes[i])
                end
                insert_face!(current_cell, 0, node_buf; is_boundary = true)
            end
            if !foundA1
                if current_cell == colA.cells[end]
                    break
                end
                for i in eachindex(colA.cells)
                    if colA.cells[i] == current_cell
                        current_cell = colA.cells[i+1]
                        ok = true
                        break
                    end
                end
            end
            it += 1
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