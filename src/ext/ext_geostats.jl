"""
    generate_perm_poro(dims, box_lengths = (1.0, 1.0, 1.0); kwarg...)

Generate porosity and permeability realizations on a box-shaped GeoStats grid.

Implementation is provided by `JutulDarcyGeoStatsExt` when `GeoStats` is loaded.
"""
function generate_perm_poro end

"""
    map_to_domain(domain, property; mapping = :mean, info_level = 0, name = missing)

Map a property sampled on a GeoStats grid or point cloud onto a JutulDarcy
reservoir domain.

Implementation is provided by `JutulDarcyGeoStatsExt` when `GeoStats` is loaded.
"""
function map_to_domain end

"""
    map_to_domain!(domain, property, name; mapping = :mean, info_level = 0, name = missing)

Mutating variant of `map_to_domain` that stores the mapped values in the named
property on the domain.
"""
function map_to_domain! end
