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
    map_realization_to_reservoir_domain(domain, realization; copy_domain = true)

Map a realization to a JutulDarcy reservoir domain.

Implementation is provided by `JutulDarcyGeoStatsExt` when `GeoStats` is loaded.
"""
function map_realization_to_reservoir_domain end
