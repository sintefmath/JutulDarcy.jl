function Jutul.triangulate_mesh(m::DataDomain{T}; kwarg...) where {T<:SimpleWell}
    return triangulate_well_domain(m; is_simple = true, kwarg...)
end

function Jutul.triangulate_mesh(m::DataDomain{T}; kwarg...) where {T<:MultiSegmentWell}
    return triangulate_well_domain(m; is_simple = false, kwarg...)
end

function triangulate_well_domain(w::DataDomain; rfactor = 10.0, outer = false, is_simple = false, flatten = true, flip = false)
    if is_simple
        error("Triangulation of SimpleWell is not implemented yet")
    else
        pts = w[:cell_centroids]
        radius = w[:radius] .* rfactor

        r = sum(radius) / length(radius)
    end
    ctri = triangulate_cylindrical_mesh(pts, radius = r, cap = :round)
    cell_index = ctri.vert_srcidx
    @info "???" cell_index size(pts)
    mapper = (
                Cells = (cell_data) -> cell_data[cell_index],
                indices = (Cells = cell_index, )
            )
    return (points = ctri.xyz', triangulation = ctri.tri', mapper = mapper)
end

