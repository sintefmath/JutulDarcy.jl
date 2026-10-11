export setup_phreeqc_model!, phreeqc_component_names, phreeqc_initial_concentrations
export close_phreeqc!, phreeqc_component_molar_masses

"""
    setup_phreeqc_model!(model, input; equilibrium_phases = String[], kwarg...)

Add fully implicit equilibrium chemistry to the existing tracer transport system.
Load `JPhreeqc` to enable this package extension (Julia 1.9 or newer).

`input` is PHREEQC solution-definition text. `equilibrium_phases` lists pure
minerals allowed to precipitate or dissolve, at saturation index zero. Tracer
unknowns are **total** component inventories, in kg component equivalent per kg
aqueous carrier, including minerals. Only the aqueous split is transported.
The component order is given by [`phreeqc_component_names`](@ref). H and O
concentrations are expressed relative to a fixed pure-water reference; use
[`phreeqc_initial_concentrations`](@ref) to obtain concentrations in this basis.

Defaults use prepared task-local workers, batches of 32 cells, exact-state
caching, and aqueous solver-guess reuse. Exact conservation derivatives are
combined with batched forward differences for active mineral sensitivities.
`temperature = 298.15` is in K; a model's `Temperature` variable takes precedence.

Minerals change `FluidVolume` by default (`update_pore_volume = true`). Molar
volumes come from the database's `PHASE_VM`, or `mineral_molar_volumes` in m³/mol,
in equilibrium-phase order. `reference_mineral_volume` is the initial volume in
m³ per cell occupied by tracked minerals, already excluded from the base pore
volume; its default of zero treats the base as mineral-free pore space. Existing
rock compressibility is retained. Flow density and viscosity remain prescribed.
See `docs/src/man/advanced/phreeqc.md` for units and supported chemistry.
"""
function setup_phreeqc_model! end

"""Return the ordered PHREEQC components transported as tracers."""
function phreeqc_component_names end

"""Return PHREEQC component formula weights in kg/mol, in tracer order."""
function phreeqc_component_molar_masses end

"""
    phreeqc_initial_concentrations(model; solution = 1)

Return component mass concentrations from a PHREEQC `SOLUTION` definition, in
tracer order. Use as initial total concentrations or injector concentrations.
Initial mineral inventory must be added to the total concentrations separately.
H and O are returned as signed excess concentrations relative to a fixed
pure-water reference, not PHREEQC's changing water-species concentration.
"""
function phreeqc_initial_concentrations end

"""Release native chemistry instances. Call after the simulation has finished."""
function close_phreeqc! end
