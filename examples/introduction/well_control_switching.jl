# # Switching a producer from rate to pressure control
# <tags: Introduction, Wells>
# Schedule a producer at two different rates with the same minimum bottom-hole
# pressure (BHP). The well meets the first request and switches to pressure
# control when the second request would take it below the BHP limit.

using Jutul, JutulDarcy

bar, day, darcy = si_units(:bar, :day, :darcy)
mesh = reservoir_mesh((3, 1, 1), (300.0, 100.0, 10.0))
domain = reservoir_domain(mesh; permeability = 0.1*darcy, porosity = 0.2)
producer = setup_well(domain, [3]; name = :Producer)
system = ImmiscibleSystem((AqueousPhase(), LiquidPhase()); reference_densities = [1000.0, 800.0])
model = setup_reservoir_model(domain, system; wells = [producer])
state0 = setup_reservoir_state(model; Pressure = 150*bar, Saturations = [1.0, 0.0])

# A pressure boundary supplies water at the opposite end of the reservoir.
# Negative rates denote production. Each entry in `forces` applies to one
# report step, with rates in m³/s and pressures in Pa.
boundary = flow_boundary_condition([1], domain, 150*bar; fractional_flow = [1.0, 0.0])
requested_rates = [-10.0, -500.0] # m³/day
minimum_bhp = 140*bar
forces = map(requested_rates) do rate
    control = Dict(:Producer => ProducerControl(TotalRateTarget(rate/day)))
    limits = Dict(:Producer => (bhp = minimum_bhp,))
    return setup_reservoir_forces(model; control, limits, bc = boundary)
end
case = JutulCase(model, [day, day], forces; state0)
result = simulate_reservoir(case; info_level = -1, error_on_incomplete = true)

# Read the operating control and the actual well response at each report time.
# When the BHP limit is active, the well produces less than requested.
operating_control = result.wells[:Producer, :control]
actual_rates = result.wells[:Producer, :rate] .* day
actual_bhp = result.wells[:Producer, :bhp] ./ bar
@assert operating_control == [:rate, :bhp] # hide
@assert isapprox(actual_rates[1], requested_rates[1]; rtol = 1e-3) # hide
@assert requested_rates[2] < actual_rates[2] < 0 # hide
@assert isapprox(actual_bhp[2], minimum_bhp/bar; rtol = 1e-3) # hide
for i in eachindex(result.time)
    println("Day $(result.time[i]/day), $(operating_control[i]) control: ",
        "requested $(requested_rates[i]) m³/day, produced $(round(actual_rates[i]; digits = 1)) m³/day ",
        "at $(round(actual_bhp[i]; digits = 1)) bar")
end
