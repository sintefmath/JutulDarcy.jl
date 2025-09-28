using JutulDarcy, Jutul, Test
Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day);
nx = ny = 5
nz = 4
dims = (nx, ny, nz)
g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
domain = reservoir_domain(g, permeability = 0.1*Darcy, porosity = 0.2)

is_simple = false
Inj = setup_well(domain, [(nx, ny, 1)], name = :Injector, reference_depth = -1.0)
Inj_w = physical_representation(Inj)
@test Inj isa DataDomain
@test Inj_w isa WellDomain
@test Inj_w.name == :Injector



Prod = setup_vertical_well(domain, 1, 1, name = :Producer, simple_well = is_simple)
Prod_shifted = setup_vertical_well(domain, 1, 1, name = :Producer, simple_well = is_simple, reference_depth = -1.0)


phases = (LiquidPhase(), VaporPhase())
rhoLS = 1000.0
rhoGS = 100.0
rhoS = [rhoLS, rhoGS] .* kg/meter^3
sys = ImmiscibleSystem(phases, reference_densities = rhoS)
model, prm = setup_reservoir_model(domain, sys, wells = [Inj, Prod], extra_out = true);

@test prm[:Injector][:FluidVolume][1] ≈ 0.392699 atol = 0.001
@test prm[:Injector][:WellIndicesThermal][1] ≈ 0.392699 atol = 0.001
@test prm[:Injector][:PerforationGravityDifference][1] ≈ 71.0982 atol = 0.001
@test prm[:Injector][:WellIndices][1] ≈ 1.18321e-12 rtol = 0.001e-12

# Simple Inj
#   :FluidVolume                  => [0.392699]
#   :WellIndicesThermal           => [27.3211]
#   :PerforationGravityDifference => [71.0982]
#   :WellIndices                  => [1.18321e-12]
#   :PhaseViscosities             => [0.001; 0.001;;]
# MS prod
# Dict{Symbol, Any} with 4 entries:
#   :PerforationGravityDifference => [0.0, 0.0, 0.0, 0.0]
#   :FluidVolume                  => [0.392699, 0.392699, 0.392699, 0.392699, 0.392699]
#   :WellIndices                  => [1.18321e-12, 1.18321e-12, 1.18321e-12, 1.18321e-12]
#   :PhaseViscosities             => [0.001 0.001 … 0.001 0.001; 0.001 0.001 … 0.001 0.001]        
