"""
    JutulDarcy.SPE10.setup_compositional_case(ncomponents = 3; kwarg...)

Set up a synthetic compositional SPE10 variant with 3 to 10 total components.
With `water = true` (the default), one component is immiscible water; the
remaining components are hydrocarbons selected
from methane through n-decane, using a Peng–Robinson equation of state with
zero binary interaction coefficients. These are benchmark fluids rather than
fluid descriptions from the original SPE10 case.

The default reservoir is liquid at 373.15 K and the SPE10 datum pressure. Its
hydrocarbon composition is 10 mol% methane, 80 mol% decane, and 10 mol%
distributed equally among the intermediate components. With only two
hydrocarbons, the composition is 10 mol% methane and 90 mol% decane. Override these
with `composition` (hydrocarbon mole fractions only) and `temperature`, or pass
an `equil` region or callback for custom initialization. Fluid properties and
relative permeabilities use the compositional model defaults. With water
enabled, the wells inject water. With `water = false`, all components are
hydrocarbons and the wells inject pure methane gas. The total number of
components remains `ncomponents` in both modes.

Other keyword arguments are forwarded to [`setup_case`](@ref), including
`layers`, `tiling`, `num_steps`, `t_total`, and `forces_arg`.

```julia
case = JutulDarcy.SPE10.setup_compositional_case(5; layers = 1:1)
result = simulate_reservoir(case)
```
"""
function setup_compositional_case(ncomponents::Integer = 3;
        water::Bool = true,
        temperature = 373.15,
        composition = missing,
        equil = missing,
        kwarg...
    )
    3 <= ncomponents <= 10 || throw(ArgumentError("ncomponents must be between 3 and 10"))
    isfinite(temperature) && temperature > 0 || throw(ArgumentError("temperature must be positive and finite"))
    hydrocarbons = ["Methane", "Ethane", "n-Propane", "n-Butane", "n-Pentane",
        "n-Hexane", "n-Heptane", "n-Octane", "n-Nonane", "n-Decane"]
    nhc = ncomponents - Int(water)
    indices = round.(Int, range(1, length(hydrocarbons); length = nhc))
    mixture = MultiComponentMixture(hydrocarbons[indices])
    eos = GenericCubicEOS(mixture, PengRobinson())
    if water
        phases = (AqueousPhase(), LiquidPhase(), VaporPhase())
        densities = [1025.0, 700.0, 0.7]
    else
        phases = (LiquidPhase(), VaporPhase())
        densities = [700.0, 0.7]
    end
    system = MultiPhaseCompositionalSystemLV(eos, phases;
        reference_densities = densities)
    if ismissing(composition)
        if nhc == 2
            composition = [0.1, 0.9]
        else
            composition = [0.1; fill(0.1/(nhc - 2), nhc - 2); 0.8]
        end
    end
    length(composition) == nhc || throw(ArgumentError("composition must have $nhc hydrocarbon mole fractions"))
    all(x -> isfinite(x) && x >= 0, composition) && sum(composition) ≈ 1.0 ||
        throw(ArgumentError("composition must contain finite, nonnegative mole fractions summing to one"))
    if ismissing(equil)
        equil = function (model)
            depths = reservoir_domain(model)[:cell_centroids][3, :]
            return EquilibriumRegion(model, 6000*si_unit(:psi), 12000*si_unit(:feet);
                temperature = temperature, composition = composition,
                woc = maximum(depths) + 100*si_unit(:feet),
                goc = minimum(depths) - 100*si_unit(:feet))
        end
    end
    return setup_case(; system = system, equil = equil, kwarg...)
end
