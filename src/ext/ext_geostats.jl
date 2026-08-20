"""
    generate_perm_poro(dims, box_lengths = (1.0, 1.0, 1.0); kwargs...)

Generate correlated porosity and permeability realizations on a box-shaped grid.

The default GeoStats-backed implementation samples a latent Gaussian process,
translates and scales it to the requested porosity mean/std, clips to the
porosity bounds, and derives permeability through `perm_from_poro`. The public
stub is defined here so the method is available from `JutulDarcy` even when the
GeoStats extension is not loaded.

Keyword arguments include:
- `nrealizations`, `seed`, `box_origin`
- `porosity_mean`, `porosity_std`, `porosity_bounds`
- `permeability_bounds`, `perm_from_poro`, and an optional custom
  `porosity_process`
"""
function generate_perm_poro end
