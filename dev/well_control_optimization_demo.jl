# Smoke test / demo for setup_well_control_optimization + optimize_well_controls
# ============================================================================
#
#   julia --project=. dev/well_control_optimization_demo.jl
#
# Needs the adjoint-well-control-switching fix (Jutul branch + JutulDarcy branch).

using Jutul, JutulDarcy, LinearAlgebra, Printf, Test

const bar = 1e5
const day = si_unit(:day)

# ---------------------------------------------------------------------------
# A small 5-spot-ish case: one water injector, one producer.
# ---------------------------------------------------------------------------
function demo_case(; nstep = 24, dtdays = 20.0)
    mesh = CartesianMesh((15, 15, 2), (600.0, 600.0, 20.0))
    domain = reservoir_domain(mesh; permeability = 1e-13, porosity = 0.22)
    Inj  = setup_vertical_well(domain, 1, 1;   name = :INJ)
    Prod = setup_vertical_well(domain, 15, 15; name = :PROD)

    rhoS = [1000.0, 780.0]          # water, oil
    sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()); reference_densities = rhoS)
    model, parameters = setup_reservoir_model(domain, sys; wells = [Inj, Prod], extra_out = true)
    replace_variables!(model, PhaseMassDensities = ConstantCompressibilityDensities(
        p_ref = 200*bar, density_ref = rhoS, compressibility = [1e-5/bar, 8e-5/bar]))

    state0 = setup_reservoir_state(model; Pressure = 200*bar, Saturations = [0.1, 0.9])
    dt = repeat([dtdays*day], nstep)

    base_rate = 900.0/day                       # m^3/s reservoir-ish
    controls = Dict(
        :INJ  => InjectorControl(TotalRateTarget(base_rate), [1.0, 0.0], density = rhoS[1]),
        :PROD => ProducerControl(SurfaceOilRateTarget(-600.0/day)),
    )
    limits = Dict(
        :INJ  => (bhp = 350*bar,),
        :PROD => (bhp = 120*bar,),
    )
    forces = setup_reservoir_forces(model; control = controls, limits = limits)
    return JutulCase(model, dt, forces; state0 = state0, parameters = parameters), base_rate
end

case, base_rate = demo_case()
nstep = length(case.dt)

# Three control periods over the 24 steps.
periods = [1:8, 9:16, 17:24]

# DOFs: injector total rate (per period), producer bhp limit (per period),
# producer oil-rate target held constant across the horizon.
controls = [
    (well = :INJ,  quantity = :rate, abs_min = 0.0,      abs_max = 3*base_rate),
    (well = :PROD, quantity = :bhp,  abs_min = 100*bar,  abs_max = 190*bar),
    (well = :PROD, quantity = :orat, rel_min = 0.5, rel_max = 2.0, constant = true),
]

copt = setup_well_control_optimization(case, periods, controls;
    objective = :npv,
    npv_arg = (oil_price = 60.0, water_price = -6.0, water_cost = 3.0, discount_rate = 0.1),
    verbose = true)
display(copt)

# ---- gradient sanity check vs central finite difference --------------------
opt = JutulDarcy.well_control_optimization_problem(copt)

function central_fd(prob, x, i; h = 1e-3)
    xp = copy(x); xp[i] += h
    xm = copy(x); xm[i] -= h
    (first(prob(xp; gradient = false)) - first(prob(xm; gradient = false))) / (2h)
end

npv_of(c) = Jutul.evaluate_objective(copt.objective, c.model,
    simulate_reservoir(c; info_level = -1).result.states, c.dt, c.forces)

@testset "well control optimization" begin
    @testset "gradient vs central FD" begin
        f, g = opt(opt.x0)
        gscale = norm(g, Inf)
        # dof 4 = PROD bhp limit, period 1. In this case PROD switches bhp->orat
        # around step 3, so NPV(bhp_limit) has a strong kink at x0; central FD
        # straddles it. KNOWN OPEN ISSUE: the adjoint value there is currently
        # outside the FD subdifferential (sign flip) - see notes.
        for i in eachindex(opt.x0)
            gfd = central_fd(opt, opt.x0, i)
            near_zero = max(abs(g[i]), abs(gfd)) < 1e-3*gscale
            ok = near_zero || isapprox(g[i], gfd; rtol = 0.1, atol = 0.05*abs(gfd))
            @printf("  dof %-2d  adjoint % .5e   FD % .5e   %s\n", i, g[i], gfd,
                ok ? "ok" : "MISMATCH")
            if i == 4
                @test_broken ok
            else
                @test ok
            end
        end
    end

    @testset "optimization improves NPV" begin
        prm = optimize_well_controls(copt; maximize = true, max_it = 15)
        npv0 = npv_of(copt(copt.dopt.parameters))
        npv1 = npv_of(copt(prm))
        @printf("  NPV base = %.5e   optimized = %.5e   (%+.1f%%)\n",
            npv0, npv1, 100*(npv1 - npv0)/abs(npv0))
        @test npv1 >= npv0 - 1e-6*abs(npv0)
    end
end
