
# An example based on the article by Zeng et al., “Effect of surfactant partitioning between gaseous phase and aqueous phase on CO₂ foam transport for 
# enhanced oil recovery”. The model runs successfully, but convergence is relatively slow. 
#

using Jutul, JutulDarcy
using HYPRE
using GLMakie
nx = 100


pth=jutul_output_path("D:/jutul_output\\case_corse_22")

Darcy, bar, kg, meter, Kelvin, day = si_units(:darcy, :bar, :kilogram, :meter, :Kelvin, :day)
ft=0.3048

cart_dims = (nx, 1, 1)
physical_dims = (1*ft, 0.01*ft, 0.01*ft)

cart_mesh = CartesianMesh(cart_dims, physical_dims)
mesh = UnstructuredMesh(cart_mesh, z_is_depth = true)

nc = number_of_cells(mesh)
perm = reshape(repeat([1.0, 1.0, 1.0], nc), 3, nc)*Darcy


poro=fill(0.25,nc)

domain_base = reservoir_domain(mesh)

mesh_presentation=Jutul.UnstructuredMesh(domain_base.representation)
domain= reservoir_domain(mesh_presentation,permeability = perm, porosity = poro, temperature = convert_to_si(120, :Celsius))

domain.data[:permeability]=(perm,Cells())

#= T0=repeat([353.15],nc)*Kelvin
domain.data[:temperature]=(T0,domain.data[:porosity][2]) =#


Injector = setup_well(domain, 1,radius = 0.0005, name = :Injector, simple_well = true,)
Producer = setup_vertical_well(domain, nc ,1 ,radius = 0.0005, name = :Producer)



physics = :kvalue

model= setup_reservoir_model(domain, :co2brine,
    wells = [Injector,Producer],
    extra_out = false,
    co2_physics = physics
);


parameters = setup_parameters(model)
parameters[:Reservoir][:PhaseViscosities] = repeat([0.045, 0.24], 1, nc)


import JutulDarcy: table_to_relperm, add_relperm_parameters!, brooks_corey_relperm
s = range(0, 1, 50)
krog_t = brooks_corey_relperm.(s, n = 4.2, kr_max=0.2,residual = 0.05,residual_total=0.1)

krog = PhaseRelativePermeability(s, krog_t, label = :og)

krg_t = brooks_corey_relperm.(s, n = 1.3, residual = 0.05, kr_max=0.94, residual_total=0.1)

krg= PhaseRelativePermeability(s, krg_t, label = :g)


#= fig, ax, plt = lines(1 .- s,  krg_t, label = "krg ")
lines!(ax, s, krog_t, label = "kro")
axislegend()
fig =#

import JutulDarcy: ReservoirRelativePermeabilities

relperm = ReservoirRelativePermeabilities(g = krg, og = krog)
replace_variables!(model, RelativePermeabilities = relperm)
add_relperm_parameters!(model);



rmodel = reservoir_model(model)

phases = JutulDarcy.get_phases(:co2brine)


#include("..//tracer_surfactant_new.jl")

import JutulDarcy: Tracers
diss=1.0
co2dens=model.models[:Reservoir].secondary_variables[:PhaseMassDensities].tab(275.79*1e5, 120+273.15)[2]
waterdens=model.models[:Reservoir].secondary_variables[:PhaseMassDensities].tab(275.79*1e5, 120+273.15)[1]

diss*=waterdens/co2dens #The form of the partition coefficient differs from that in Zeng et al., so a conversion is needed.


Tracer1= Tracers.SurfactantTracer(rmodel.system,k_d=diss,Concentration_limit="Liquid",C_max=0.1,)#D=(gas=1e-5,liquid=1e-11))
#Tracer2=MultiPhaseTracer(rmodel.system)


add_tracers_to_model!(model,Tracer1)

Tracers.set_surfactant_model!(model)


fm=Tracers.FmFactor(epdry=500.0, fmdry=0.25, fmsurf=2.0, epsurf=1.0, fommb=500.0)

replace_variables!(rmodel, FmFactor = fm)



pv = pore_volume(model, parameters)

inj_rate=5.66336*1e-5/86400

total_time=0.5*sum(pv)/inj_rate/86400 #day



dt = [total_time*day]#CO2 and surf

inj_target= TotalRateTarget(inj_rate)#m3/s
I_ctrl = InjectorControl(inj_target, [0, 1.0], density = co2dens,tracers=[2.5/co2dens])


controls= Dict()
controls[:Injector] = I_ctrl


#bc = flow_boundary_condition(nc, domain, 275.79*bar,convert_to_si(120, :Celsius), fractional_flow = [1.0,0.0])

bhp_target = BottomHolePressureTarget(275.79*bar)
P_ctrl= ProducerControl(bhp_target)

controls[:Producer] = P_ctrl

forces= setup_reservoir_forces(model, control =controls)   #bc=bc)



push!(rmodel.output_variables, :PhaseViscosities)
push!(rmodel.output_variables, :FmFactor)
push!(rmodel.output_variables, :ClVolumetricConcentration)
push!(rmodel.output_variables, :RelativePermeabilities)


parameters = setup_parameters(model)
state0 = setup_reservoir_state(model, Pressure = 275.79*bar, Saturations = [1.0, 0.0],OverallMoleFractions = [1.0, 0.0])
ws, states = simulate_reservoir(state0, model, dt, 
    parameters = parameters,
    forces=forces,
    max_timestep = 2e-3/86400*day,
    info_level=2,
    restart = false,
    output_path=pth,
);

