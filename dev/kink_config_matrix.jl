# Config matrix for the zero-gradient / kink puzzle
# ==================================================
#
# From dev/switching_kink_investigation.jl we know, for the single-DOF setup
# (PROD bhp limit, period 1 only, scaler = missing):
#   * NPV(b1) is piecewise smooth with slope ~ -0.096 per Pa around 120 bar and
#     large negative slopes at 160-185 bar,
#   * there is an isolated upward "reset" of NPV exactly at b1 == 120 bar (the
#     value of the periods-2/3 limits), suggesting a force-equality artifact
#     (e.g. timestep selection reacting to "forces changed" between periods),
#   * the adjoint gradient is IDENTICALLY ZERO everywhere - even at 170 bar
#     where all 8 period-1 report steps operate on the bhp limit.
#
# The 3-period demo config gave a correct nonzero adjoint (-0.0969/Pa). This
# script isolates which ingredient kills the gradient, and verifies the
# force-equality/timestep hypothesis for the reset at 120.
#
# Run:  julia --project=. dev/kink_config_matrix.jl

using Jutul, JutulDarcy, Printf, LinearAlgebra

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
npv_arg = (oil_price = 60.0, water_price = -6.0, water_cost = 3.0, discount_rate = 0.1)

function make_prob(spec; backend_arg = missing)
    copt = setup_well_control_optimization(case, periods, [spec];
        objective = :npv, npv_arg = npv_arg, verbose = false)
    kw = ismissing(backend_arg) ? NamedTuple() : (backend_arg = backend_arg,)
    prob = JutulDarcy.well_control_optimization_problem(copt; kw...)
    return copt, prob
end

# x-vector for a physical magnitude vector `vals` (one per DOF entry),
# accounting for the linear_limits scaler if present.
function x_of(spec, vals)
    if get(spec, :scaler, missing) === :linear_limits
        lo = spec.abs_min; hi = spec.abs_max
        return [(v - lo)/(hi - lo) for v in vals]
    else
        return collect(Float64, vals)
    end
end

# convert dNPV/dx (component i) to physical dNPV/db in 1/Pa
function to_phys(spec, g)
    if get(spec, :scaler, missing) === :linear_limits
        return g / (spec.abs_max - spec.abs_min)
    else
        return g
    end
end

function probe(name, spec, vals; icomp = 1, h_phys = 0.2*bar, backend_arg = missing)
    copt, prob = make_prob(spec; backend_arg = backend_arg)
    x = x_of(spec, vals)
    f, g = prob(x)
    g_phys = to_phys(spec, g[icomp])
    # one-sided FD in physical units (both sides), consistent objective via prob
    dx = x_of(spec, [h_phys])[1] - x_of(spec, [0.0])[1]  # h in x-units
    xp = copy(x); xp[icomp] += dx
    xm = copy(x); xm[icomp] -= dx
    fp = first(prob(xp; gradient = false))
    fm = first(prob(xm; gradient = false))
    fwd = to_phys(spec, (fp - f)/dx)
    bwd = to_phys(spec, (f - fm)/dx)
    @printf("%-42s  f = %.6e\n", name, f)
    @printf("    adjoint = %+.5e   FD fwd = %+.5e   FD bwd = %+.5e   (per Pa)\n",
        g_phys, fwd, bwd)
    flush(stdout)
    return (f = f, g = g_phys, fwd = fwd, bwd = bwd)
end

b_eval = 170.0*bar   # all 8 period-1 report steps on the bhp limit here
println("=" ^ 78)
println("A) Which configuration kills the adjoint?  (PROD bhp limit at 170 bar)")
println("   true slope from sweep: ~ -1.8e+00 per Pa (NPV drops ~4.4e5 per 2.5 bar)")
println("=" ^ 78)

probe("1. periods=[1], scaler=missing",
    (well = :PROD, quantity = :bhp, abs_min = b_lo, abs_max = b_hi,
     periods = [1], scaler = missing), [b_eval])

probe("2. periods=[1], scaler=:linear_limits",
    (well = :PROD, quantity = :bhp, abs_min = b_lo, abs_max = b_hi,
     periods = [1], scaler = :linear_limits), [b_eval])

probe("3. periods=[1,2,3], scaler=missing (d/d period1)",
    (well = :PROD, quantity = :bhp, abs_min = b_lo, abs_max = b_hi,
     scaler = missing), [b_eval, 120*bar, 120*bar])

probe("4. periods=[1,2,3], scaler=:linear_limits",
    (well = :PROD, quantity = :bhp, abs_min = b_lo, abs_max = b_hi,
     scaler = :linear_limits), [b_eval, 120*bar, 120*bar])

probe("5. constant=true, scaler=missing (all periods)",
    (well = :PROD, quantity = :bhp, abs_min = b_lo, abs_max = b_hi,
     constant = true, scaler = missing), [b_eval])

println()
println("=" ^ 78)
println("B) Backend variants for config 1 (the zero-gradient config)")
println("=" ^ 78)

probe("1a. dense (use_sparsity=false, di_sparse=false)",
    (well = :PROD, quantity = :bhp, abs_min = b_lo, abs_max = b_hi,
     periods = [1], scaler = missing), [b_eval];
    backend_arg = (use_sparsity = false, di_sparse = false, single_step_sparsity = false, do_prep = true))

probe("1b. deps_ad = :di",
    (well = :PROD, quantity = :bhp, abs_min = b_lo, abs_max = b_hi,
     periods = [1], scaler = missing), [b_eval];
    backend_arg = (deps_ad = :di,))

println()
println("=" ^ 78)
println("C) Sanity: a TARGET dof with periods=[1] (PROD orat, period 1)")
println("=" ^ 78)
orat0 = 600.0/day
probe("6. PROD orat target, periods=[1], scaler=missing",
    (well = :PROD, quantity = :orat, abs_min = 0.5*orat0, abs_max = 2*orat0,
     periods = [1], scaler = missing), [orat0]; h_phys = 0.01*orat0)

println()
println("=" ^ 78)
println("D) The reset at b1 == 120 bar: ministep-count hypothesis")
println("=" ^ 78)
copt0, _ = make_prob((well = :PROD, quantity = :bhp, abs_min = b_lo, abs_max = b_hi,
    periods = [1], scaler = missing))
for b in [119.8, 119.9, 119.95, 120.0, 120.05, 120.1, 120.2] .* bar
    prm = deepcopy(copt0.dopt.parameters)
    prm["PROD"]["bhp"] = [b]
    c = copt0.setup_function(prm, missing)
    states, reports = simulate(c; output_substates = true, info_level = -1)
    n_mini = [length(filter(m -> m[:success], r[:ministeps])) for r in reports]
    G = (model, state, dt, step_info, forces) -> JutulDarcy.npv_objective(
        model, state, dt, step_info, forces;
        injectors = [:INJ], producers = [:PROD], timesteps = c.dt, npv_arg...)
    f = Jutul.evaluate_objective(G, c.model, states, c.dt, c.forces)
    @printf("  b1 = %8.3f bar   NPV = %.7e   substates = %3d   ministeps/report step = %s\n",
        b/bar, f, length(states), join(n_mini, ","))
    flush(stdout)
end

println("DONE")
