function Jutul.triangulate_mesh(m::DataDomain{T}; kwarg...) where {T<:SimpleWell}
    return triangulate_well_domain(m; is_simple = true, kwarg...)
end

function Jutul.triangulate_mesh(m::DataDomain{T}; kwarg...) where {T<:MultiSegmentWell}
    return triangulate_well_domain(m; is_simple = false, kwarg...)
end

function triangulate_well_domain(w::DataDomain; rfactor = 10.0, outer = false, is_simple = false, flatten = true, flip = false)
    if is_simple
        pts = w[:perforation_centroids]
    else
        pts = w[:cell_centroids]
    end
    # if size(pts, 2) == 1
    #     pts = repeat(pts, 1, 2)
    #     pt_remap = [1, 1]
    # else
    #     pt_remap = 1:size(pts, 2)
    # end
    pt_remap = collect(1:size(pts, 2))
    if size(pts, 1) == 3
        pushfirst!(pt_remap, 1)
        pts = hcat(pts[:, 1], pts)
        pts[3, 1] -= 10.0
    end

    radius = w[:radius] .* rfactor
    r = sum(radius) / length(radius)

    ctri = triangulate_cylindrical_mesh(pts, radius = r, cap = :round)
    cell_index = pt_remap[ctri.vert_srcidx]
    if is_simple
        cell_index .= 1
    end
    mapper = (
        Cells = (cell_data) -> cell_data[cell_index],
        indices = (Cells = cell_index, )
    )
    return (points = ctri.xyz', triangulation = ctri.tri', mapper = mapper)
end

