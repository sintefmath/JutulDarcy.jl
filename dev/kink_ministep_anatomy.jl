# E5/E6: ministep anatomy across the 120-bar discontinuity + dense-mode test
# ===========================================================================
#
# E5: The one-sided jump of NPV at b1 = 120 bar (period-1 bhp floor crossing
#     the periods-2/3 value) must come from a discrete change in the simulation
#     path. Dump the ministep structure and per-substate operating controls on
#     both sides, using the SAME simulator configuration as the optimization
#     problem (setup_reservoir_simulator with output_substates).
#
# E6: Frozen-sparsity test: does dense differentiation (use_sparsity = false,
#     di_sparse = false) avoid the "first gradient at inactive limit -> all
#     later gradients zero" failure?
#
# Run:  julia --project=. dev/kink_ministep_anatomy.jl

using Jutul, JutulDarcy, Printf

const bar = 1e5
const day = si_unit(:day)

function demo_case(; nstep = 24, dtdays = 20.0)
    mesh = CartesianMesh((15, 15, 2), (600.0, 600.0, 20.0))
    domain = reservoir_domain(mesh; permeability = 1e-13, porosity = 0.22)
    Inj  = setup_vertical_well(domain, 1, 1;   name = :INJ)
    Prod = setup_vertical_well(domain, 15, 15; name = :PROD)
    rhoS = [1000.0, 780.0]
    sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()); reference_densities = rhoS)
    model, parameters = setup_reservoir_model(domain, sys; wells = [Inj, Prod], extra_out = true)
    replace_variables!(model, PhaseMassDensities = ConstantCompressibilityDensities(
        p_ref = 200*bar, density_ref = rhoS, compressibility = [1e-5/bar, 8e-5/bar]))
    state0 = setup_reservoir_state(model; Pressure = 200*bar, Saturations = [0.1, 0.9])
    dt = repeat([dtdays*day], nstep)
    controls = Dict(
        :INJ  => InjectorControl(TotalRateTarget(900.0/day), [1.0, 0.0], density = rhoS[1]),
        :PROD => ProducerControl(SurfaceOilRateTarget(-600.0/day)),
    )
    limits = Dict(:INJ => (bhp = 350*bar,), :PROD => (bhp = 120*bar,))
    forces = setup_reservoir_forces(model; control = controls, limits = limits)
    return JutulCase(model, dt, forces; state0 = state0, parameters = parameters)
end

case = demo_case()
periods = [1:8, 9:16, 17:24]
b_lo, b_hi = 95*bar, 185*bar
spec = (well = :PROD, quantity = :bhp, abs_min = b_lo, abs_max = b_hi,
    periods = [1], scaler = missing)
npv_arg = (oil_price = 60.0, water_price = -6.0, water_cost = 3.0, discount_rate = 0.1)

copt = setup_well_control_optimization(case, periods, [spec];
    objective = :npv, npv_arg = npv_arg, verbose = false)

function anatomy(b)
    prm = deepcopy(copt.dopt.parameters)
    prm["PROD"]["bhp"] = [b]
    c = copt.setup_function(prm, missing)
    # Same solver setup as the optimization problem uses
    sim, cfg = setup_reservoir_simulator(c;
        output_substates = true, info_level = -1, end_report = false)
    res = simulate!(sim, c.dt; state0 = c.state0, parameters = c.parameters,
        forces = c.forces, config = cfg)
    states, reports = res.states, res.reports
    nsub = length(states)
    # ministep dt structure per report step
    mini = [[m[:dt] for m in r[:ministeps] if m[:success]] for r in reports]
    # substate operating controls + bhp for PROD
    ops = Symbol[]
    bhps = Float64[]
    for s in states
        cfgw = s[:Facility][:WellGroupConfiguration]
        ctrl = cfgw.operating_controls[:PROD]
        push!(ops, ctrl isa DisabledControl ? :disabled :
            JutulDarcy.translate_target_to_symbol(ctrl.target))
        push!(bhps, s[:Facility][:BottomHolePressure][
            JutulDarcy.get_well_position(sim.model[:Facility].domain, :PROD)])
    end
    return (nsub = nsub, mini = mini, ops = ops, bhps = bhps)
end

println("="^78)
println("E5) Ministep anatomy across the b1 = 120 bar discontinuity")
println("="^78)
for b in [119.90, 119.95, 119.99, 120.00, 120.01, 120.05, 121.00] .* bar
    a = anatomy(b)
    @printf("\n-- b1 = %8.3f bar   substates = %d\n", b/bar, a.nsub)
    for (i, m) in enumerate(a.mini)
        if length(m) > 1 || i <= 2 || 8 <= i <= 10
            @printf("   step %2d: %d ministeps, dt(days) = %s\n", i, length(m),
                join([@sprintf("%.3f", x/day) for x in m], ", "))
        end
    end
    println("   PROD op ctrl per substate: ", join(first.(string.(a.ops)), ""))
    println("   PROD bhp per substate (bar): ",
        join([@sprintf("%.2f", x/bar) for x in a.bhps[1:min(8, end)]], ", "), " ...")
    flush(stdout)
end

println()
println("="^78)
println("E6) Dense-mode frozen-pattern test (gradient at 100 bar FIRST, then 170)")
println("    Sparse reference: second gradient came back 0.0 (wrong, true -1.817)")
println("="^78)
prob_dense = JutulDarcy.well_control_optimization_problem(copt;
    backend_arg = (use_sparsity = false, di_sparse = false, single_step_sparsity = false, do_prep = true))
f1, g1 = prob_dense([100.0*bar])
@printf("  dense grad @100 = %+.5e (expect ~0)\n", g1[1])
f2, g2 = prob_dense([170.0*bar])
@printf("  dense grad @170 = %+.5e (expect -1.81710e+00 if dense avoids the freeze)\n", g2[1])

println("DONE")
