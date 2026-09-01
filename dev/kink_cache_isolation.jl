# Isolating the two caching artifacts found by kink_config_matrix.jl
# ===================================================================
#
# Established so far (dev/kink_config_matrix.jl):
#   * A FRESH problem object gives a correct adjoint (-1.817 vs central FD
#     -1.816 at 170 bar) for every DOF configuration.
#   * The first script (dev/switching_kink_investigation.jl) reused ONE problem
#     object: its first gradient eval was at 100 bar (limit inactive, gradient
#     legitimately zero) - and every subsequent gradient came back 0.0, even at
#     170 bar. -> suspicion: adjoint storage / sparsity frozen at first eval.
#   * Fresh simulations show NPV(b1) is FLAT near 120 bar, while the cached-
#     simulator sweep showed a spurious smooth slope (-0.096/Pa) and a sawtooth
#     reset at 120.00 -> suspicion: forward-evaluation hysteresis via state
#     (e.g. WellGroupConfiguration / timestep selector) carried across
#     evaluations in the cached simulator.
#
# This script isolates both.
#
# Run:  julia --project=. dev/kink_cache_isolation.jl

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

function fresh_prob()
    copt = setup_well_control_optimization(case, periods, [spec];
        objective = :npv,
        npv_arg = (oil_price = 60.0, water_price = -6.0, water_cost = 3.0, discount_rate = 0.1),
        verbose = false)
    return JutulDarcy.well_control_optimization_problem(copt)
end

# ---------------------------------------------------------------------------
println("="^78)
println("E1) Forward-evaluation hysteresis in the cached simulator")
println("    Same prob object; f(120 bar) evaluated after different predecessors.")
println("="^78)
let prob = fresh_prob()
    seq = [120.0, 185.0, 120.0, 100.0, 120.0, 120.0, 170.0, 120.0] .* bar
    prev = "-"
    for b in seq
        f = first(prob([b]; gradient = false))
        @printf("  f(%5.1f bar) = %.7e    (previous eval: %s)\n", b/bar, f, prev)
        prev = @sprintf("%.1f bar", b/bar)
        flush(stdout)
    end
end

# ---------------------------------------------------------------------------
println()
println("="^78)
println("E2) Stale adjoint storage: gradient correctness vs evaluation order")
println("    Reference (fresh prob, gradient at 170 first): adjoint = -1.817e0/Pa")
println("="^78)
let prob = fresh_prob()
    println("  -- prob A: gradient at 100 bar FIRST (limit inactive), then 170:")
    f1, g1 = prob([100.0*bar])
    @printf("    grad @100 = %+.5e (expect ~0, correct)\n", g1[1])
    f2, g2 = prob([170.0*bar])
    @printf("    grad @170 = %+.5e (fresh-prob reference: -1.81710e+00)\n", g2[1])
    flush(stdout)
end
let prob = fresh_prob()
    println("  -- prob B: gradient at 170 bar FIRST, then 100, then 170 again:")
    f1, g1 = prob([170.0*bar])
    @printf("    grad @170 = %+.5e (expect -1.81710e+00)\n", g1[1])
    f2, g2 = prob([100.0*bar])
    @printf("    grad @100 = %+.5e (expect ~0)\n", g2[1])
    f3, g3 = prob([170.0*bar])
    @printf("    grad @170 = %+.5e (repeat - stale storage would corrupt this)\n", g3[1])
    flush(stdout)
end

# ---------------------------------------------------------------------------
println()
println("="^78)
println("E3) The sawtooth: does the reset location depend on sweep direction?")
println("    Ascending vs descending fine sweep with ONE cached prob each.")
println("="^78)
bs = collect(range(119.6*bar, 120.4*bar; length = 17))
let prob = fresh_prob()
    f_up = [first(prob([b]; gradient = false)) for b in bs]
    prob2 = fresh_prob()
    f_dn = [first(prob2([b]; gradient = false)) for b in reverse(bs)]
    reverse!(f_dn)
    @printf("  %-10s %-16s %-16s %-12s\n", "b (bar)", "f ascending", "f descending", "diff")
    for (b, fu, fd) in zip(bs, f_up, f_dn)
        @printf("  %8.3f  %.7e   %.7e   %+9.1f\n", b/bar, fu, fd, fu - fd)
    end
    flush(stdout)
end

# ---------------------------------------------------------------------------
println()
println("="^78)
println("E4) True local shape near 120 bar from FRESH problems (no cache reuse)")
println("="^78)
for b in [118.0, 119.0, 119.5, 120.0, 120.5, 121.0, 122.0] .* bar
    prob = fresh_prob()
    f = first(prob([b]; gradient = false))
    @printf("  f(%7.2f bar) = %.7e   (fresh prob+simulator)\n", b/bar, f)
    flush(stdout)
end

println("DONE")
