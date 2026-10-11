# Reactive transport with PHREEQC

Loading [JPhreeqc.jl](https://github.com/simulkade/JPhreeqc.jl) enables
`JutulDarcyJPhreeqcExt` on Julia 1.9 or newer. The extension adds equilibrium
speciation and pure-mineral precipitation/dissolution to the tracer equations.
Its native library currently supports x86-64 Windows and Linux.

```julia
using Jutul, JutulDarcy, JPhreeqc

mesh = CartesianMesh((20, 1), (20.0, 1.0))
domain = reservoir_domain(mesh, porosity = 0.2, permeability = 1e-12)
model = setup_reservoir_model(domain, SinglePhaseSystem(AqueousPhase()))
input = """
SOLUTION 1
    temp 25
    pH 9
    units mol/kgw
    Na 0.01
    Cl 0.01
    Ca 0.002
    C(4) 0.004
END
"""
setup_phreeqc_model!(model, input, equilibrium_phases = ["Calcite"])
names = phreeqc_component_names(model)
c0 = phreeqc_initial_concentrations(model)
state0 = setup_reservoir_state(model, Pressure = 101325.0,
    TracerConcentrations = c0)

try
    result = simulate_reservoir(state0, model, [1.0], info_level = -1)
finally
    close_phreeqc!(model)
end
```

## Conserved totals and the phase split

Transport **components** (elements and charge), rather than individual aqueous
species. PHREEQC resolves species and pure minerals locally. Its component
ordering is authoritative; always obtain it with `phreeqc_component_names`.

`TracerConcentrations` contains total component equivalents, including mineral
inventory, in kg/kg aqueous carrier. `PhreeqcAqueousConcentrations` is the mobile
split in the same units. `PhreeqcImmobileConcentrations` is their difference.
`TracerMasses` contains total equivalents in kg. The carrier mass is
`FluidVolume * Saturations[aqueous_phase] * PhaseMassDensities[aqueous_phase]`.
Thus the conservation law accumulates the **total** inventory while face and
well-perforation fluxes transport the **aqueous** concentration.

H and O use a fixed pure-water reference, rather than subtracting the changing
PHREEQC water species. With component molar masses `M` in kg/mol, the reference
mass fractions are `2M_H/(2M_H + M_O)` and `M_O/(2M_H + M_O)`.
The H/O tracer values are total elemental mass fractions minus those constants.
This avoids transporting large water-background values while preserving a
consistent elemental basis across speciation and precipitation. H, O and Charge
tracers can be negative. Charge uses PHREEQC's formula-weight convention.
`phreeqc_initial_concentrations` performs this conversion automatically.

To initialize existing minerals, add their elemental inventories to `c0`.
For example, for `m` mol of calcite per kg aqueous carrier, add
`m*M_Ca`, `m*M_C` and `3m*M_O` to the respective totals, using
`phreeqc_component_molar_masses(model)`. Do not add a nonzero mineral reservoir
to the PHREEQC input. The equilibrium phase list supplies zero-amount seeds;
all available mass comes from conserved totals. Precipitated material can
dissolve on later steps as the transported totals and equilibrium change.

Injection concentrations come from a solution definition, for example
`phreeqc_initial_concentrations(model, solution = 0)`. Pass this vector as
`InjectorControl(..., tracers = cin)`. Simple wells carry aqueous concentrations;
reservoir production samples the aqueous split rather than the total inventory.

## Mineral volumes and pore space

`PhreeqcMineralConcentrations` contains mineral amounts in mol/kg aqueous
carrier, in `equilibrium_phases` order. `PhreeqcMineralVolume` contains their
combined specific solid volume in m³/kg carrier. The extension retrieves actual
phase amounts from PHREEQC, so phases with identical elemental formulas remain
distinct. Elemental aqueous totals are reconstructed from mineral stoichiometry
and total inventories, preserving component balances exactly.

Minerals change `FluidVolume` by default. Their molar volumes come from
[`PHASE_VM`](https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-61.htm)
in the database, converted from cm³/mol to m³/mol. Override them with
`mineral_molar_volumes = [...]`, in equilibrium-phase order and **m³/mol**.
Missing positive volumes produce an error when pore-volume feedback is enabled.
Set `update_pore_volume = false` for fixed pore space; databases without mineral
volumes then report zero specific volume unless volumes are supplied.

Let `V_b` denote the base pore volume supplied by the rock model,
`V_s,ref` the volume of tracked minerals already excluded from the specified
base pore volume, and `ν_s` their current specific volume. The local closure is

```math
V = V_b + V_{s,ref} - \rho_a S_a V\nu_s
  = \frac{V_b + V_{s,ref}}{1 + \rho_a S_a\nu_s}.
```

This expression uses the current aqueous density and saturation and carries all
derivatives into flow and tracer accumulation. It requires no lagged porosity or
timestep history. For constant base volume, the existing parameter becomes
`StaticFluidVolume`. An existing pressure-dependent `FluidVolume` secondary
variable becomes `PhreeqcUnreactedFluidVolume`; rock compressibility is retained.
Configure rock compressibility **before** adding chemistry.

`reference_mineral_volume` accepts a scalar or a vector of **m³ per cell**,
defaulting to zero. Zero treats the domain pore volume as available space before
tracked minerals are added. If the specified porosity already excludes existing
tracked minerals, supply their occupied volume. For prescribed initial amounts
`n_p` in mol/kg carrier, this is `V_0*ρ_a,0*S_a,0*sum(n_p .* molar_volumes)`.
Initial elemental inventories must still be included in tracer totals. The
reference volume can also be supplied through the
`PhreeqcReferenceMineralVolume` parameter when setting up simulation parameters.
Effective porosity is `FluidVolume / bulk_cell_volume`; permeability and flow
PVT properties remain prescribed.

## Derivatives and repeatability

`PhreeqcPhaseSplit` is the chemistry secondary variable; aqueous concentrations,
mineral amounts and specific volume are extracted from its common output.
Its dependencies include total tracers, pressure, and the model's temperature if present.
Temperature is in K and pressure is absolute Pa; the extension converts them
to PHREEQC's Celsius and atmospheres. Native cells use a normalized 1 L volume
and 1 kg/L carrier density with kg/kg solution units. Pure-mineral equilibrium
is homogeneous in extensive quantities, so its concentrations and specific
mineral volumes are independent of carrier density, saturation and pore volume.
Flow density and saturation enter the pore-volume closure and accumulation
directly. PHREEQC-calculated solution volume and density do not replace flow PVT.

PHREEQC builds a chemical Newton Jacobian internally, but its public
[C interface](https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqcrm-html/RM__interface__C_8h.html)
and JPhreeqc do not expose equilibrium sensitivities. Species concentrations and
activity coefficients alone do not provide those derivatives.

The extension uses exact identity derivatives for spectator-element balances
and for cells without active minerals in the interior of an undersaturated
region. It differentiates only mineral amounts numerically and reconstructs
aqueous derivatives from elemental conservation. Density derivatives of the
normalized phase split are exactly zero; constant temperature is not perturbed.
Cross-component, pressure and variable-temperature sensitivities propagate
through Jutul AD. The immobile split, mineral volume, reactive pore volume and
accumulation retain their derivatives. Forward differences are the economical
default; `difference = :central` provides a validation mode, using forward
differences at nonnegative inventory bounds.
Derivatives are local to the active mineral set; they are one-sided at phase
appearance/disappearance. Higher derivatives are not supported.

Every evaluation and perturbation resets mineral inventories to zero and supplies
the complete conserved totals. Aqueous activities are retained as nonlinear
initial guesses (`reuse_speciation = true`), which reduces solver work without
retaining physical history. Set `reuse_speciation = false` to restore the full
baseline on every flash. Recomputed equilibria agree to solver accuracy; exact
cache hits are bit-identical. Failed native solves force a full restore on the
next call. Newton iterations, timestep cuts and restarts cannot advance mineral
history. Native errors return to Julia with cell-range context rather
than aborting the process. Saved primary totals reconstruct the equilibrium
split on restart; native handles are not restart data.

## Defaults and threading

Start Julia with multiple threads, for example `julia --threads=8 --project`.
The extension preconfigures up to `Threads.nthreads()` native instances, each
with one native thread. A task checks out a prepared worker once and binds it in
scoped `task_local_storage`; nested property calls retrieve that instance and its
reusable buffers. The binding ends before the worker returns to the pool. Task
migration and Jutul's own threaded secondary-variable loop therefore do not
share a native handle concurrently. No databases or instances are created in
the property-update loop.
Native calls operate on batches of 32 cells with reusable numeric buffers.
Smaller models use fewer workers. `threads` and `batch_size` are configurable.

Chemistry file/screen output and individual-species storage are disabled.
Native `PRINT` controls also disable input echo, status updates and `USER_GRAPH`;
these are not all disabled by turning off chemistry printing alone.
Minimal in-memory selected output saves only mineral amounts and saturation
indices. Native load balancing is disabled for the single-thread instances.
PHREEQC's internal `KNOBS -numerical_derivatives` option is false by default;
its native analytic derivatives are used where supported by the activity model.
This controls PHREEQC's own nonlinear solve, not the external phase-split
sensitivities, and does not expose its internal Jacobian through the C API.
Native convergence tolerance is `1e-12` to resolve derivative perturbations.
`relative_step = 1e-5` controls finite-difference steps. Exact-state caches reuse
results and derivatives only when every chemistry input matches; there is no
approximate chemistry skipping. Cache memory scales with the number of cells
times the square of the component count.

`INCREMENTAL_REACTIONS true` accelerates multi-step native batch kinetics by
reusing the previous interval's result. It changes the meaning of listed
reaction/time increments and supplies no material benefit to this coupling's
single-step, zero-time equilibrium `RunCells` calls. It is therefore not forced
as a performance default. Native kinetic options such as `-cvode` and
`-step_divide` likewise do not apply to the current equilibrium map.
See the [USGS reaction-step documentation](https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-19.htm).

`concentration_scale = 1e-5` gives reactive tracers stricter convergence than
passive tracers. With the standard nonlinear tolerance of `1e-3`, the tracer
criterion bounds the component residual over a timestep by approximately
`1e-8` kg equivalent per m³ pore volume. Change the scale or solver tolerances
deliberately for the chemistry and concentrations of interest.

## Supported scope

The current coupling is fully implicit equilibrium chemistry in a present
aqueous phase, with pure minerals at saturation index zero. Kinetics, exchange,
surface sites, gas phases, solid solutions and irreversible mineral constraints
need additional history variables and are rejected in solution input. Mineral
inventories modify effective pore volume; flow density, viscosity, permeability
and thermal mineral properties are not recalculated from chemistry. The tracer
formulation assumes a prescribed carrier fluid with dilute transported components.
The same mineral list applies to all reservoir cells; solutions and initial
inventories can vary by cell. Existing passive tracer sets, multisegment wells,
GPU chemistry, MPI and NLDD subdomain chemistry are not supported. Set up
chemistry before creating simulator storage; close it after all simulations
using the model have finished. Model copies share the synchronized runtime.

## Validation

`test/phreeqc.jl` checks an independent IPhreeqc calcite reference, elemental
stoichiometry, serial/threaded agreement, repeated evaluations, all phase-split
partials, heterogeneous native batches and padding, aqueous multiphase
upwinding, flow and total-inventory conservation in a closed flow simulation,
task-local worker reuse and cleanup, exact inactive/spectator derivatives,
precipitation/dissolution volume changes, reactive volume partials and rock
compressibility, Jutul's threaded secondary-variable loop, and error recovery.
It requires JPhreeqc in the active test environment:

```julia
using JutulDarcy, JPhreeqc
include(joinpath(dirname(dirname(pathof(JutulDarcy))), "test", "phreeqc.jl"))
```

`examples/validation/phreeqc_column.jl` demonstrates mineral dissolution in a
flowing column with injection and production wells. Run it with Julia threads
enabled. `examples/validation/phreeqc_batch.jl` initializes existing calcite
inventory and compares equilibrium at two temperatures. The full test suite
treats JPhreeqc validation as optional when its
platform-specific dependency is absent.

The unit and component conventions follow the
[PhreeqcRM C interface documentation](https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqcrm-html/RM__interface__C_8h.html).
