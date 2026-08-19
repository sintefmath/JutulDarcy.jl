"""
    setup_perm_poro_realizations(dims, box_lengths; kwarg...)

Generate porosity/permeability realizations over a box.

Implementation is provided by `JutulDarcyGeoStatsExt` when `GeoStats` is loaded.
"""
function generate_perm_poro end

"""
    map_to_domain(domain, property; box, mapping = :mean, copy_domain = false, info_level = 0)

Map a box-defined property vector onto a JutulDarcy reservoir domain.

Implementation is provided by `JutulDarcyGeoStatsExt` when `GeoStats` is loaded.
"""
function map_to_domain end

"""
    map_to_domain!(domain, property, name; box, mapping = :mean, copy_domain = false, info_level = 0)

Mutating variant of `map_to_domain` that stores the mapped values in the named
property on the domain.
"""
function map_to_domain! end
