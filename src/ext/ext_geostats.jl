"""
    realization = generate_perm_poro(mesh; kwarg...)
    realizations = generate_perm_poro(mesh; nrealizations = 10, kwarg...)
    realization = generate_perm_poro(domain::DataDomain; kwarg...)
    realization = generate_perm_poro(geostats_domain; kwarg...)

Generate spatially correlated porosity and permeability realizations using
[GeoStats.jl](https://juliaearth.github.io/GeoStats.jl). Requires that GeoStats
is loaded (`using GeoStats`).

A standard Gaussian field `z` is sampled from a Gaussian process (by default
with a spherical covariance), transformed to porosity and clamped to
`porosity_bounds`, and permeability is computed from porosity with
`perm_from_poro`. Each realization is returned as a `Dict` with entries
`:porosity` and `:permeability`, and also `:geotable` (the GeoStats realization
on the domain where it was sampled) if `extra_out = true`. A single `Dict` is
returned if `nrealizations = 1`, otherwise a vector of `Dict`s.

# Where the field is sampled
For a Jutul mesh (or a `DataDomain` with a mesh), `method` decides where the
Gaussian field is sampled:
- `:ijk`: On the logical IJK grid of the mesh, using the average cell sizes as
  grid spacing (see `ijk_spacing`). Each cell gets the value of its logical
  cell, so the correlation follows the layering of e.g. corner-point meshes,
  and no geometric mapping is needed. Default for meshes with IJK structure.
- `:centroids`: At the cell centroids. Works for any mesh. Default for meshes
  without IJK structure.
- `:box`: On a uniform `CartesianGrid` covering the bounding box of the mesh
  with `box_dims` cells, after which porosity and permeability are mapped to
  the mesh with [`map_property_to_cells`](@ref) (arithmetic averaging for
  porosity, `perm_mapping` for permeability).

A GeoStats domain (e.g. a `CartesianGrid` or a `GeometrySet`) can also be given
directly, in which case the values are ordered as the elements of that domain.

# Keyword arguments
- `nrealizations = 1`: Number of realizations.
- `rng = Random.default_rng()`: Random number generator.
- `seed = nothing`: If given, use `Random.Xoshiro(seed)` as generator.
- `porosity_mean = 0.2`, `porosity_std = 0.05`: Porosity is
  `porosity_mean + porosity_std*z`.
- `porosity_transform = nothing`: Optional function `z -> porosity` replacing
  the linear transform above, e.g. a normal-score transform to a target
  distribution.
- `porosity_bounds = (0.05, 0.95)`: Porosity is clamped to these bounds.
- `correlation_range`: Correlation range in all directions. Defaults to a
  third of the horizontal and vertical extent of the domain.
- `horizontal_correlation_range`, `vertical_correlation_range`: Override the
  horizontal and vertical ranges separately.
- `porosity_process = nothing`: Custom GeoStats field process for `z`, or a
  function `(rng, n) -> z` returning `n` values. Overrides the ranges above.
- `simulation_method = nothing`: GeoStats simulation method (e.g. `FFTSIM()`,
  `LUSIM()` or `SEQSIM()`). By default GeoStats picks one based on the domain
  and the correlation range. `FFTSIM()` is much faster for large grids, but is
  only available for `:ijk` and `:box` and should only be used when the
  correlation range is small compared to the grid.
- `perm_from_poro = kozeny_carman_permeability`: Function
  `(porosity; bounds) -> permeability`.
- `permeability_bounds = (1e-20, Inf)`: Bounds passed to `perm_from_poro`.
- `kv_kh = nothing`: If given, permeability is returned as a 3 by n matrix
  with the vertical permeability scaled by this ratio.
- `method`, `ijk_spacing`, `box_dims`, `perm_mapping = :geometric`: See
  above and [`map_property_to_cells`](@ref).
- `info_level = 0`: Print a summary if positive.
- `extra_out = false`: Include the GeoStats realization as `:geotable`.
"""
function generate_perm_poro end

"""
    realization = generate_perm_poro!(domain::DataDomain; kwarg...)

Generate a single porosity and permeability realization with
[`generate_perm_poro`](@ref) and store it in `domain` as `:porosity` and
`:permeability`. Returns the realization `Dict`. Requires that GeoStats is
loaded (`using GeoStats`).
"""
function generate_perm_poro! end
