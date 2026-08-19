"""
    kozeny_carman_permeability(porosity; kwarg...)

Compute permeability from porosity.

Implementation is provided by `JutulDarcyGeoStatsExt` when `GeoStats` is loaded.
"""
function kozeny_carman_permeability end

"""
    setup_perm_poro_realizations(dims, box_lengths; kwarg...)

Generate porosity/permeability realizations over a box.

Implementation is provided by `JutulDarcyGeoStatsExt` when `GeoStats` is loaded.
"""
function setup_perm_poro_realizations end

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

"""
    map_realization_to_reservoir_domain(domain, realization; copy_domain = true)

Backward-compatible wrapper around `map_to_domain` for legacy realization-based inputs.
"""
function map_realization_to_reservoir_domain end
