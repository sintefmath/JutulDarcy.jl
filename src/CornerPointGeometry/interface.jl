
"""
    set_cpgrid_geometry!(domain)
    set_cpgrid_geometry!(domain, grid_section)

Replace the interior computed geometry `domain` with a variant computed from the
corner nodes of the GRID section of a corner-point file. If you do not pass a
`grid_section` this function will only change the geometry if `:cpgrid_corners`
has first been added with the internal routine
[`add_cpgrid_corners_to_reservoir!`](@ref).

Subsequent use of geometry to compute transmissibilities will be consistent with
typical reservoir simulators as these use corner-point geometry for their
internal calculations instead of the "true" geometry of the resulting polyogonal
cells.
"""
function set_cpgrid_geometry!(domain, grid_section = missing)
    geometry = cpgrid_geometry(domain, grid_section)
    if !isnothing(geometry)
        domain[:cell_centroids] .= geometry.cell_centroids
        domain[:face_centroids] .= geometry.face_centroids
        domain[:areas] .= geometry.face_areas
        domain[:normals] .= geometry.face_normals
    end
    return geometry
end
