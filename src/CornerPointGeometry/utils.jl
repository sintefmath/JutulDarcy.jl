using Jutul.StaticArrays: SVector

"""
    cpgrid_cell_corners(grid, ijk)

Reconstruct the eight logical corners of deck cell `ijk` from the `COORD` and
`ZCORN` keywords in a parsed `GRID` section.
"""
function cpgrid_cell_corners(grid::AbstractDict, ijk, dims = Tuple(grid["cartDims"]), coord = reshape(grid["COORD"], 6, :)')
    nx, ny, _ = dims
    zcorn = grid["ZCORN"]
    T = promote_type(eltype(grid["COORD"]), eltype(zcorn))
    linear = GeoEnergyIO.CornerPointGrid.ijk_to_linear(ijk..., dims)
    corners = ntuple(8) do ix
        i = (ix - 1) % 2
        j = ((ix - 1) ÷ 2) % 2
        k = ((ix - 1) ÷ 4) % 2
        zix = GeoEnergyIO.CornerPointGrid.corner_index(linear, (i, j, k), dims)
        line = GeoEnergyIO.CornerPointGrid.get_line(coord, ijk[1] + i, ijk[2] + j, nx + 1, ny + 1)
        return SVector{3, T}(GeoEnergyIO.CornerPointGrid.interp_coord(line[1], line[2], zcorn[zix]))
    end
    return SVector{8, SVector{3, T}}(corners)
end

@inline cpgrid_corner_index(i, j, k) = 1 + i + 2j + 4k

function cpgrid_face_corner_indices(direction, side)
    if direction == 1
        return SVector(
            cpgrid_corner_index(side, 0, 0), cpgrid_corner_index(side, 1, 0),
            cpgrid_corner_index(side, 1, 1), cpgrid_corner_index(side, 0, 1)
        )
    elseif direction == 2
        return SVector(
            cpgrid_corner_index(0, side, 0), cpgrid_corner_index(0, side, 1),
            cpgrid_corner_index(1, side, 1), cpgrid_corner_index(1, side, 0)
        )
    elseif direction == 3
        return SVector(
            cpgrid_corner_index(0, 0, side), cpgrid_corner_index(1, 0, side),
            cpgrid_corner_index(1, 1, side), cpgrid_corner_index(0, 1, side)
        )
    end
    throw(ArgumentError("Expected logical face direction 1, 2, or 3, got $direction."))
end

function cpgrid_area_normal(points)
    center = sum(points)/length(points)
    area_normal = zero(first(points))
    previous = points[end]
    for point in points
        area_normal += cross(previous - center, point - center)/2
        previous = point
    end
    return area_normal
end

"""
    add_cpgrid_corners_to_reservoir!(domain, grid)

Cache the eight static CPGRID corners for every active cell as the
`:cpgrid_corners` cell property of `domain`. The parsed `GRID` section must
contain `COORD`, `ZCORN`, and `cartDims`.
"""
function add_cpgrid_corners_to_reservoir!(domain::DataDomain, grid::AbstractDict)
    for key in ("COORD", "ZCORN", "cartDims")
        haskey(grid, key) || throw(ArgumentError("GRID section is missing required keyword $key."))
    end
    mesh = Jutul.physical_representation(domain)
    nc = Jutul.number_of_cells(mesh)
    T = promote_type(eltype(grid["COORD"]), eltype(grid["ZCORN"]))
    dims = Tuple(grid["cartDims"])
    coord = reshape(grid["COORD"], 6, :)'
    corners = Vector{SVector{8, SVector{3, T}}}(undef, nc)
    for cell in 1:nc
        corners[cell] = cpgrid_cell_corners(grid, Jutul.cell_ijk(mesh, cell), dims, coord)
    end
    domain[:cpgrid_corners, Jutul.Cells()] = corners
    return domain
end

"""
    cpgrid_geometry(domain, grid_section = missing)
    cpgrid_geometry(domain)

Build the Flow-style shared geometry from the cached `:cpgrid_corners` cell
property. Returns `nothing` when CPGRID corners have not been added to `domain`.
You can pass the `GRID` section as the second argument to automatically add the
corners if they are not already present. Only ZCORN/COORD formats are supported.
"""
function cpgrid_geometry(domain, grid_section = missing)
    if !ismissing(grid_section)
        add_cpgrid_corners_to_reservoir!(domain, grid_section)
    end
    haskey(domain, :cpgrid_corners) || return nothing
    mesh = Jutul.physical_representation(domain)
    nc = Jutul.number_of_cells(mesh)
    nf = Jutul.number_of_faces(mesh)
    neighbors = domain[:neighbors]
    size(neighbors, 2) == nf || throw(ArgumentError("Mesh has $nf faces but domain neighborship has $(size(neighbors, 2))."))

    corners = domain[:cpgrid_corners]
    T = eltype(eltype(first(corners)))
    cell_centroids = zeros(T, 3, nc)
    for cell in 1:nc
        cell_centroids[:, cell] .= sum(corners[cell])/8
    end
    face_centroids = zeros(T, 3, nf)
    face_areas = zeros(T, nf)
    face_normals = zeros(T, 3, nf)

    face_dir = JutulDarcy.get_ijk_face_dir(mesh, neighbors)

    for face in 1:nf
        left, right = neighbors[:, face]
        direction = face_dir[face]
        left_ijk = Jutul.cell_ijk(mesh, left)
        right_ijk = Jutul.cell_ijk(mesh, right)
        logical_sign = sign(right_ijk[direction] - left_ijk[direction])
        logical_sign != 0 || throw(ArgumentError("Face $face has no I/J/K change between cells $left and $right."))
        side = logical_sign > 0 ? 1 : 0
        points = corners[left][cpgrid_face_corner_indices(direction, side)]
        # Use half face distances from the cell-side CPGRID centers
        area_normal = cpgrid_area_normal(mesh.node_points[mesh.faces.faces_to_nodes[face]])
        area = norm(area_normal)
        face_centroids[:, face] .= sum(points)/4
        face_areas[face] = area
        if area < 1e-12
            # Pinched faces carry no flux. Keep a finite placeholder normal.
            face_normals[direction, face] = logical_sign
        else
            face_normals[:, face] .= area_normal/area
        end
    end

    return (
        cell_centroids = cell_centroids,
        face_centroids = face_centroids,
        face_areas = face_areas,
        face_normals = face_normals
    )
end

