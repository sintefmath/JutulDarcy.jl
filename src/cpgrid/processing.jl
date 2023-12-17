
@enum CPGRID_PILLAR_AB_INTERSECTION begin
    AB_RANGES_MATCH # A and B have the same top and bottom points on this pillar
    AB_OVERLAP_A_FIRST # A and B overlap, with A on top and then B starting after A but continuing further along
    AB_OVERLAP_B_FIRST
    A_CONTAINS_B # A fully contains B and is larger in both upwards and downwards direction
    B_CONTAINS_A
    TOP_MATCHES_A_LONG # The top point matches and A goes further down than B
    TOP_MATCHES_B_LONG
    BOTTOM_MATCHES_A_LONG # The bottom point matches and A goes further up than B
    BOTTOM_MATCHES_B_LONG
    DISTINCT_A_ABOVE # The ranges do not overlap and A is above B
    DISTINCT_B_ABOVE
end


function find_cell_bounds(cell, line)
    get_pos(i) = line.cellpos[i]:(line.cellpos[i+1]-1)
    get_cells(i) = @view line.cells[get_pos(i)]

    n = length(line.cellpos)-1
    start = 0
    for i in 1:n
        if cell in get_cells(i)
            start = i
            break
        end
    end
    @assert start > 0 "$cell start point not found in line $line"
    stop = 0
    for i in (n:-1:start)
        if cell in get_cells(i)
            stop = i
            break
        end
    end
    @assert stop > 0 "$cell end point not found in line $line"
    return (start, stop)
end

function cell_top_bottom(cells, line1, line2; check = true)
    # out = Vector{Tuple{Tuple{Int, Int}, Tuple{Int, Int}}}()
    T = @NamedTuple{cell::Int64, line1::Tuple{Int64, Int64}, line2::Tuple{Int64, Int64}}
    prev1 = prev2 = 1
    out = T[]
    for cell in cells
        pos1 = find_cell_bounds(cell, line1)
        pos2 = find_cell_bounds(cell, line2)
        line1_skip = prev1 != pos1[1]
        line2_skip = prev2 != pos2[1]
        if line1_skip || line2_skip
            # Add boundary ghost cell
            ghost1 = (prev1, pos1[1])
            ghost2 = (prev2, pos2[1])
            push!(out, (cell = 0, line1 = ghost1, line2 = ghost2))
        end
        prev1 = pos1[end]
        prev2 = pos2[end]

        # Positions in line1, positions in line 2
        push!(out, (cell = cell, line1 = pos1, line2 = pos2))
    end
    if check
        start1 = start2 = 1
        for (i, o) in enumerate(out)
            @assert o.line1[1] == start1
            @assert o.line2[1] == start2

            start1 = o.line1[2]
            start2 = o.line2[2]
        end
        @assert start1 == out[end].line1[2]
        @assert start2 == out[end].line2[2]
    end
    return out
end

function add_vertical_cells_from_overlaps!(extra_node_lookup, F, nodes, cell_pairs, overlaps, l1, l2)
    complicated_overlap_categories = (DISTINCT_A_ABOVE, DISTINCT_B_ABOVE)
    node_pos = Int[]
    function global_node_index(l, local_node)
        return l.nodes[local_node]
    end
    function global_node_point(l, local_node::Int)
        return nodes[global_node_index(l, local_node)]
    end

    function global_node_point_and_index(l, local_node)
        ix = global_node_index(l, local_node)
        return (nodes[ix], ix)
    end

    sizehint!(node_pos, 6)
    for (overlap, cell_pair) in zip(overlaps, cell_pairs)
        cell_a, cell_b = cell_pair

        # @info "Starting" overlap cell_pair
        cat1 = overlap.line1.category
        cat2 = overlap.line2.category

        # TODO: Figure out sign
        edge1 = overlap.line1.overlap
        edge2 = overlap.line2.overlap
        n1 = length(edge1)
        n2 = length(edge2)

        @assert n1 == 0 || maximum(edge1) <= length(l1.nodes)
        @assert n2 == 0 || maximum(edge2) <= length(l2.nodes)

        line1_distinct = cat1 == DISTINCT_A_ABOVE || cat1 == DISTINCT_B_ABOVE
        line2_distinct = cat2 == DISTINCT_A_ABOVE || cat2 == DISTINCT_B_ABOVE

        empty!(node_pos)
        if n1 + n2 < 3 && false
            # @warn "Skipping" n1 n2 cat1 cat2
            continue
        elseif cat1 == cat2 == AB_RANGES_MATCH  && false
            # TODO: This can probably include more categories.
            # @info "Simple matching!"
            # Need line1 lookup here.
            handle_matching_overlap!(node_pos, edge1, edge2, l1, l2)
        elseif true
            handle_generic_intersections!(node_pos, extra_node_lookup, nodes, cell_a, cell_b, l1, l2, overlap, global_node_point)
    
        elseif line1_distinct && line2_distinct
            handle_no_overlapping_edges(node_pos, extra_node_lookup, nodes, cell_a, cell_b, l1, l2, overlap, global_node_point)
        elseif line1_distinct
            @assert n1 == 0
            a_range = overlap.line1.range_a
            b_range = overlap.line1.range_b
            @warn "Handle1" a_range b_range cat1

            handle_one_side_distinct!(node_pos, extra_node_lookup, nodes, cat1, cell_a, cell_b, l1, l2, edge2, a_range, b_range, global_node_point)
        elseif line2_distinct
            @assert n2 == 0 "n2 = $n2 for $edge2"
            a_range = overlap.line2.range_a
            b_range = overlap.line2.range_b

            @warn "Handle2" a_range b_range cat2
            # TODO: Check that this part is actually ok.
            # handle_one_side_distinct!(node_pos, extra_node_lookup, nodes, cat2, cell_a, cell_b, l2, l1, edge1, a_range, b_range, global_node_point_and_index)
            handle_one_side_distinct!(node_pos, extra_node_lookup, nodes, cat2, cell_a, cell_b, l2, l1, edge1, a_range, b_range, global_node_point)
        else
            # Some degree of mismatch
            handle_matching_overlap_with_intersections!(node_pos, extra_node_lookup, nodes, cell_a, cell_b, l1, l2, overlap, global_node_point)
        end
        if cell_a == cell_b
            # If both cells are the same, we are dealing with a mirrored
            # column pair at the boundary.
            cell_a = 0
        end
        if length(node_pos) > 2
            F(cell_a, cell_b, node_pos)
        end
    end
    return nothing
end

function handle_matching_overlap!(node_pos, edge1, edge2, l1, l2)
    for local_node_ix in edge1
        push!(node_pos, l1.nodes[local_node_ix])
    end
    for local_node_ix in reverse(edge2)
        push!(node_pos, l2.nodes[local_node_ix])
    end
end

function handle_matching_overlap_with_intersections!(node_pos, extra_node_lookup, nodes, cell_a, cell_b, l1, l2, overlap, global_node_point)
    # Full ranges for line 1
    edge1 = overlap.line1.overlap
    range1_a = overlap.line1.range_a
    range1_b = overlap.line1.range_b
    # Full ranges for line 2
    edge2 = overlap.line2.overlap
    range2_a = overlap.line2.range_a
    range2_b = overlap.line2.range_b

    # Fist traverse down left vertical edge, then find potential crossing point
    # at lower horizontal edge, then up right vertical edge and finally
    # potential crossing point at top horizontal edge.

    # 1. Left vertical edge overlap
    for local_node_ix in edge1
        push!(node_pos, l1.nodes[local_node_ix])
    end

    # 2. Potential crossing point at lower horizontal edge
    a1_l = range1_a[end]
    b1_l = range1_b[end]
    a2_l = range2_a[end]
    b2_l = range2_b[end]
    a_above_left_b_above_right = a1_l < b1_l && a2_l > b2_l
    b_above_left_a_above_right = a1_l > b1_l && a2_l < b2_l
    if a_above_left_b_above_right || b_above_left_a_above_right
        # Crossing between two top lines
        l1_p1 = global_node_point(l1, a1_l)
        l1_p2 = global_node_point(l2, a2_l)
        l2_p1 = global_node_point(l1, b1_l)
        l2_p2 = global_node_point(l2, b2_l)
        bottom_edge_node = handle_crossing_node!(extra_node_lookup, nodes, l1_p1, l1_p2, l2_p1, l2_p2)
        push!(node_pos, bottom_edge_node)
    end

    # 3. Reverse traversal of right edge
    for local_node_ix in reverse(overlap.line2.overlap)
        push!(node_pos, l2.nodes[local_node_ix])
    end

    # 4. Finally potential crossing point at top horizontal edge
    a1_t = range1_a[1]
    b1_t = range1_b[1]
    a2_t = range2_a[1]
    b2_t = range2_b[1]

    a_above_left_b_above_right = a1_t < b1_t && a2_t > b2_t
    b_above_left_a_above_right = a1_t > b1_t && a2_t < b2_t

    if a_above_left_b_above_right || b_above_left_a_above_right
        l1_p1 = global_node_point(l1, a1_t)
        l1_p2 = global_node_point(l2, a2_t)
        l2_p1 = global_node_point(l1, b1_t)
        l2_p2 = global_node_point(l2, b2_t)
        bottom_edge_node = handle_crossing_node!(extra_node_lookup, nodes, l1_p1, l1_p2, l2_p1, l2_p2)
        push!(node_pos, bottom_edge_node)
    end
end

function handle_no_overlapping_edges(node_pos, extra_node_lookup, nodes, cell_a, cell_b, l1, l2, overlap, global_node_point)
    # Full ranges for line 1
    edge1 = overlap.line1.overlap
    range1_a = overlap.line1.range_a
    range1_b = overlap.line1.range_b
    # Full ranges for line 2
    edge2 = overlap.line2.overlap
    range2_a = overlap.line2.range_a
    range2_b = overlap.line2.range_b

    first_and_last(x) = (first(x), last(x))
    a1_t, a1_l = first_and_last(range1_a)
    a2_t, a2_l = first_and_last(range2_a)
    b1_t, b1_l = first_and_last(range1_b)
    b2_t, b2_l = first_and_last(range2_b)

    # Lower nodes
    l1_p1_l = global_node_point(l1, a1_l)
    l1_p2_l = global_node_point(l2, a2_l)
    l2_p1_l = global_node_point(l1, b1_l)
    l2_p2_l = global_node_point(l2, b2_l)
    # Top nodes
    l1_p1_t = global_node_point(l1, a1_t)
    l1_p2_t = global_node_point(l2, a2_t)
    l2_p1_t = global_node_point(l1, b1_t)
    l2_p2_t = global_node_point(l2, b2_t)

    # lower to lower
    n_ll = handle_crossing_node!(extra_node_lookup, nodes, l1_p1_l, l1_p2_l, l2_p1_l, l2_p2_l)
    push!(node_pos, n_ll)

    # lower [1] to top [2]
    n_lt = handle_crossing_node!(extra_node_lookup, nodes, l1_p1_l, l1_p2_l, l2_p1_t, l2_p2_t)
    push!(node_pos, n_lt)

    # top to top 
    n_tt = handle_crossing_node!(extra_node_lookup, nodes, l1_p1_t, l1_p2_t, l2_p1_t, l2_p2_t)
    push!(node_pos, n_tt)

    # top [1] to lower [2]
    n_tl = handle_crossing_node!(extra_node_lookup, nodes, l1_p1_t, l1_p2_t, l2_p1_l, l2_p2_l)
    push!(node_pos, n_tl)
end

function handle_generic_intersections!(node_pos, extra_node_lookup, nodes, cell_a, cell_b, l1, l2, overlap, global_node_point)
    # Full ranges for line 1
    edge1 = overlap.line1.overlap
    range1_a = overlap.line1.range_a
    range1_b = overlap.line1.range_b
    # Full ranges for line 2
    edge2 = overlap.line2.overlap
    range2_a = overlap.line2.range_a
    range2_b = overlap.line2.range_b

    first_and_last(x) = (first(x), last(x))
    a1_t, a1_l = first_and_last(range1_a)
    a2_t, a2_l = first_and_last(range2_a)
    b1_t, b1_l = first_and_last(range1_b)
    b2_t, b2_l = first_and_last(range2_b)

    function pos_diff(a, b)
        return (a < b, a == b)
    end

    # Top-top matching
    at_before_bt_1, matching_tt_1 = pos_diff(a1_t, b1_t)
    at_before_bt_2, matching_tt_2 = pos_diff(a2_t, b2_t)
    # Low-low matching
    al_before_bl_1, matching_ll_1 = pos_diff(a1_l, b1_l)
    al_before_bl_2, matching_ll_2 = pos_diff(a2_l, b2_l)

    # Low-high matching (a > b)
    # at_before_bl_1, matching_tl_1 = pos_diff(a1_t, b1_l)
    # at_before_bl_2, matching_tl_2 = pos_diff(a2_t, b2_l)

    # High-low matching (a < b)
    # al_before_bt_1, matching_tl_1 = pos_diff(a1_t, b1_l)
    # al_before_bt_2, matching_tl_2 = pos_diff(a2_t, b2_l)

    # Top nodes
    l1_p1_t = global_node_point(l1, a1_t)
    l1_p2_t = global_node_point(l2, a2_t)
    l2_p1_t = global_node_point(l1, b1_t)
    l2_p2_t = global_node_point(l2, b2_t)

    line_top_a = (l1_p1_t, l1_p2_t)
    line_top_b = (l2_p1_t, l2_p2_t)

    # Lower nodes
    l1_p1_l = global_node_point(l1, a1_l)
    l1_p2_l = global_node_point(l2, a2_l)
    l2_p1_l = global_node_point(l1, b1_l)
    l2_p2_l = global_node_point(l2, b2_l)

    line_low_a = (l1_p1_l, l1_p2_l)
    line_low_b = (l2_p1_l, l2_p2_l)

    # 2 edge reversed first.
    #TODO: Decide order here
    # Four conditionals, each potentially adding a point
    # 1_top crossing 2_top (reversal a/b top over pair)
    if at_before_bt_1 != at_before_bt_2 && !(matching_tt_1 || matching_tt_2)
        @info "Adding top"
        @assert a1_t != b1_t
        @assert a2_t != b2_t

        n_tt = handle_crossing_node!(extra_node_lookup, nodes, line_top_a, line_top_b)
        push!(node_pos, n_tt)
    else
        @info "Skipping top"
    end

    if length(edge1) == 0
        # Check which one is the lowest
        # 1_top crossing 2_low (reversal ab_top_1 vs ab_low_1 over pair)
        if at_before_bt_1
            # A is above B at left edge. This extra node is formed from A lower
            # cutting B upper.
            n_tl = handle_crossing_node!(extra_node_lookup, nodes, line_low_a, line_top_b)
        else
            n_tl = handle_crossing_node!(extra_node_lookup, nodes, line_top_a, line_low_b)
        end
        push!(node_pos, n_tl)
    else
        for local_node_ix in edge1
            push!(node_pos, l1.nodes[local_node_ix])
        end
    end

    if al_before_bl_1 != al_before_bl_2 && !(matching_ll_1 || matching_ll_2)
        @info "Adding bottom"
        @assert a1_l != b1_l
        @assert a2_l != b2_l

        # 1_low crossing 2_low (reversal a/b low over pair)
        n_ll = handle_crossing_node!(extra_node_lookup, nodes, line_low_a, line_low_b)
        push!(node_pos, n_ll)
    else
        @info "Skipping bottom"
    end

    if length(edge2) == 0
        # 1_low crossing 2_top (reversal ab_low_1 vs ab_top_1 over pair)
        n_lt = handle_crossing_node!(extra_node_lookup, nodes, line_top_a, line_low_b)

        if at_before_bt_2
            # A is above B at left edge. This extra node is formed from A lower
            # cutting B upper.
            n_tl = handle_crossing_node!(extra_node_lookup, nodes, line_low_a, line_top_b)
        else
            n_tl = handle_crossing_node!(extra_node_lookup, nodes, line_top_a, line_low_b)
        end

        push!(node_pos, n_lt)
    else
        for local_node_ix in reverse(edge2)
            push!(node_pos, l2.nodes[local_node_ix])
        end
    end
    @info "Done" node_pos
end



function handle_one_side_distinct!(node_pos, extra_node_lookup, nodes, cat_self, cell_a, cell_b, l_self, l_other, other_edge, a_range, b_range, global_node_point)
    @info "Hey..."
    # l_self: Line where cells are distinct
    # l_other: Line where cells are not distinct (i.e. overlapping to some extent.)
    self_a_top = a_range[1]
    self_a_bottom = a_range[end]

    self_b_top = b_range[1]
    self_b_bottom = b_range[end]

    # Points on the other edge, limited to overlap
    if cat_self == DISTINCT_A_ABOVE
        upper_self = self_a_bottom
        lower_self = self_b_top
    else
        @assert cat_self == DISTINCT_B_ABOVE
        upper_self = self_b_bottom
        lower_self = self_a_top
    end
    @assert upper_self < lower_self "Error: $upper_self >= $lower_self"
    upper_pt_self = global_node_point(l_self, upper_self)
    lower_pt_self = global_node_point(l_self, lower_self)
    # These are the same no matter hat
    upper_pt_other = global_node_point(l_other, first(other_edge))
    lower_pt_other = global_node_point(l_other, last(other_edge))

    new_node_ix = handle_crossing_node!(extra_node_lookup, nodes, upper_pt_self, lower_pt_other, lower_pt_self, upper_pt_other)
    push!(node_pos, new_node_ix)
    for local_node_ix in other_edge
        push!(node_pos, l_other.nodes[local_node_ix])
    end
end

function handle_crossing_node!(extra_node_lookup, nodes, line_a, line_b)
    pt = find_crossing_node(line_a, line_b)
    return cpgrid_get_or_add_crossing_node!(extra_node_lookup, nodes, pt)
end

function find_crossing_node(l1, l2)
    @assert length(l1) == 2 == length(l2)
    return find_crossing_node(l1[1], l1[2], l2[1], l2[2])
end

function find_crossing_node(x1, x2, x3, x4)
    # Line (x1 -> x2) crossing line (x3 -> x4)
    a = x2 - x1
    b = x4 - x3
    c = x3 - x1
    axb = cross(a, b)
    nrm = sum(axb.^2)
    if nrm == 0.0 
        @warn "Bad intercept: $((x1, x2)) to $((x3, x4))" a b c
        error()
    end
    x_intercept = x1 + a*(dot(cross(c, b), axb)/nrm)
    return x_intercept
end

function cpgrid_get_or_add_crossing_node!(extra_node_lookup, nodes, pt)
    if haskey(extra_node_lookup, pt)
        ix = extra_node_lookup[pt]
    else
        push!(nodes, pt)
        ix = length(nodes)
        extra_node_lookup[pt] = ix
    end
    return ix
end

function traverse_column_pair(col_a, col_b, l1, l2)
    cell_pairs = Tuple{Int, Int}[]
    # overlaps = Tuple{UnitRange{Int64}, UnitRange{Int64}}[]
    # TODO: Make concrete.
    overlaps = []

    function find_end(a, b, s::Symbol)
        end_a = a[s][2]
        end_b = b[s][2]
        if end_a == end_b
            # Reached end of both.
            return (end_a, true, true)
        elseif end_a < end_b
            # Reached end of a but not b
            return (end_a, true, false)
        else
            # Reached end of b but not a
            return (end_b, false, true)
        end
    end

    ord_a = cell_top_bottom(col_a.cells, l1, l2)
    ord_b = cell_top_bottom(col_b.cells, l1, l2)
    function get_local_line(pos, is_line1::Bool)
        if is_line1
            d = pos.line1
        else
            d = pos.line2
        end
        start, stop = d
        return (start, stop)
    end

    function gen_category(t::CPGRID_PILLAR_AB_INTERSECTION, rng, range_a, range_b)
    return (
        category = t,
        overlap = rng,
        range_a = range_a,
        range_b = range_b
    )
end

    function determine_overlap(pos_a, pos_b, is_line1)
        a_start, a_stop = get_local_line(pos_a, is_line1)
        b_start, b_stop = get_local_line(pos_b, is_line1)
        return determine_cell_overlap_inside_line(a_start, a_stop, b_start, b_stop)
    end
    for pos_a in ord_a
        for pos_b in ord_b
            # @info "Cell pair: $((pos_a.cell, pos_b.cell))"

            t1, overlap_1, range_1a, range_1b = determine_overlap(pos_a, pos_b, true)
            t2, overlap_2, range_2a, range_2b = determine_overlap(pos_a, pos_b, false)

            t_equal = t1 == t2
            if t_equal && t1 in (DISTINCT_A_ABOVE, DISTINCT_B_ABOVE)
                continue
            end
            # Unless A > B for both pillars, or B > A for both pillars,
            # we have found a face.
            push!(overlaps,
                (
                line1 = gen_category(t1, overlap_1, range_1a, range_1b),
                line2 = gen_category(t2, overlap_2, range_2a, range_2b)
                )
            )
            push!(cell_pairs, (pos_a.cell, pos_b.cell))
        end
    end
    return (cell_pairs, overlaps)
end

function split_overlaps_into_interior_and_boundary(cell_pairs, overlaps)
    cell_is_boundary(x) = x < 1
    interior_pairs = similar(cell_pairs, 0)
    interior_overlaps = similar(overlaps, 0)

    boundary_pairs = similar(interior_pairs)
    boundary_overlaps = similar(interior_overlaps)

    for (cell_pair, overlap) in zip(cell_pairs, overlaps)
        l, r = cell_pair
        l_bnd = cell_is_boundary(l)
        r_bnd = cell_is_boundary(r)
        if l_bnd && r_bnd
            # Two inactive cells, can be skipped.
            continue
        elseif l_bnd || r_bnd
            push!(boundary_pairs, cell_pair)
            push!(boundary_overlaps, overlap)
        else
            push!(interior_pairs, cell_pair)
            push!(interior_overlaps, overlap)
        end
    end
    return (interior_pairs, interior_overlaps, boundary_pairs, boundary_overlaps)
    # return (interior = (interior_pairs, interior_overlaps), boundary = (boundary_pairs, boundary_overlaps))
end
##
function determine_cell_overlap_inside_line(a_start, a_stop, b_start, b_stop)
    a_range = a_start:a_stop
    b_range = b_start:b_stop
    # @info "Checking overlap" a_range b_range a_start a_stop b_start b_stop
    a_starts_before_b = a_start < b_start
    b_starts_before_a = b_start < a_start

    if a_range == b_range
        out = (AB_RANGES_MATCH, a_range)
    elseif a_starts_before_b && b_stop < a_stop
        out = (A_CONTAINS_B, b_range)
    elseif b_starts_before_a && a_stop < b_stop
        out = (B_CONTAINS_A, a_range)
    elseif a_start == b_start
        # Take the shortest of the two ranges.
        if a_stop > b_stop
            out = (TOP_MATCHES_A_LONG, b_range)
        else
            out = (TOP_MATCHES_B_LONG, a_range)
        end
    elseif a_stop == b_stop
        if a_starts_before_b
            # A starts earlier and is therefore longer
            out = (BOTTOM_MATCHES_A_LONG, b_range)
        else
            out = (BOTTOM_MATCHES_B_LONG, a_range)
        end
    elseif b_start > a_stop
        out = (DISTINCT_A_ABOVE, 0:-1)
    elseif a_start > b_stop
        out = (DISTINCT_B_ABOVE, 0:-1)
    elseif a_starts_before_b && b_stop > a_stop
        out = (AB_OVERLAP_A_FIRST, b_start:a_stop)
    elseif b_starts_before_a && a_stop > b_stop
        out = (AB_OVERLAP_B_FIRST, a_start:b_stop)
    else
        error("Programming error")
    end
    # @info "Result:" out
    return (out[1], out[2], a_range, b_range)
end