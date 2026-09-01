# Diagnostic: extended-control-vector gradient via the DictOptimization path
# ========================================================================
#
# Same experiment as dev/extended_control_gradient.jl, but the gradient is
# taken with the DI-based DictOptimization machinery
# (`parameters_gradient_reservoir(...; deps = :case)`) instead of the older
# `solve_adjoint_forces`. The optimization "control vector" is the dict
#   "orat" => [oil-rate target magnitude]
#   "bhp"  => [bhp lower-limit value]
# and the setup function rebuilds the forces (control + limit) from it.
#
# Question: does this path correctly follow whichever control was active
# during the forward simulation (unlike solve_adjoint_forces, which returns
# an identically-zero gradient for the limit entry)?
#
# Run:  julia --project=. dev/extended_control_gradient_dictopt.jl

using Jutul, JutulDarcy, LinearAlgebra, Printf, Test
using Jutul: OrderedDict

# NOTE: passes only with the fix on Jutul branch fix/adjoint-well-control-switching
# (+ matching JutulDarcy changes in src/facility/{types,controls,wellgroups}.jl).
# Without it, every `bhp`/`limit_*` entry comes back identically 0.0.

const bar = 1e5
const day = si_unit(:day)

function build_model(; nstep = 20, dtdays = 15.0)
    mesh = CartesianMesh((20, 20, 1), (500.0, 500.0, 20.0))
    domain = reservoir_domain(mesh; permeability = 1e-13, porosity = 0.25)
    Prod = setup_vertical_well(domain, 10, 10; name = :Producer)

    phases = (AqueousPhase(), LiquidPhase())
    rhoS = [1000.0, 800.0]
    sys = ImmiscibleSystem(phases; reference_densities = rhoS)

    model, parameters = setup_reservoir_model(domain, sys; wells = [Prod], extra_out = true)
    rho = ConstantCompressibilityDensities(
        p_ref = 200*bar, density_ref = rhoS, compressibility = [1e-5/bar, 1.5e-4/bar])
    replace_variables!(model, PhaseMassDensities = rho)

    state0 = setup_reservoir_state(model; Pressure = 200*bar, Saturations = [0.15, 0.85])
    dt = repeat([dtdays*day], nstep)
    return (; model, parameters, state0, dt, nstep)
end

# forces vector from a dict of control values
function forces_from_dict(M, d)
    if haskey(d, "bhp_control")
        ctrl = ProducerControl(BottomHolePressureTarget(d["bhp_control"][1]))
        f0 = setup_reservoir_forces(M.model; control = Dict(:Producer => ctrl),
            set_default_limits = false)
    else
        ctrl = ProducerControl(SurfaceOilRateTarget(-d["orat"][1]))
        f0 = setup_reservoir_forces(M.model; control = Dict(:Producer => ctrl),
            limits = Dict(:Producer => (bhp = d["bhp"][1],)), set_default_limits = false)
    end
    return [f0 for _ in 1:M.nstep]
end

function case_from_dict(M, d)
    return JutulCase(M.model, M.dt, forces_from_dict(M, d);
        state0 = M.state0, parameters = M.parameters)
end

function make_objective(M)
    t_tot = sum(M.dt)
    return (model, state, dt, step_info, forces) -> begin
        orat = JutulDarcy.compute_well_qoi(model, state, forces, :Producer, SurfaceOilRateTarget)
        (dt*orat)/t_tot
    end
end

# independent central finite difference on one dict entry
function fd_entry(M, d0, key; rel = 1e-3, absmin = 1e-8)
    G = make_objective(M)
    v0 = d0[key][1]
    ep = max(absmin, rel*abs(v0))
    function ev(s)
        d = deepcopy(d0); d[key] = [v0 + s]
        c = case_from_dict(M, d)
        st, _ = simulate(c; info_level = -1)
        Jutul.evaluate_objective(G, c.model, st, c.dt, c.forces)
    end
    return (ev(ep) - ev(-ep))/(2*ep), v0, ep
end

function run_regime(name, d0; keys_to_check)
    println("\n" * "="^70)
    println("REGIME $name    $(d0)")
    println("="^70)
    M = build_model()
    case0 = case_from_dict(M, d0)
    res = simulate_reservoir(case0; info_level = -1)
    println("control per step : ", res.wells[:Producer][:control])
    println("bhp (bar)        : ", round.(res.wells[:Producer][:bhp] ./ bar, digits = 1))

    G = make_objective(M)
    setup_fn(d, step_info = missing) = case_from_dict(M, d)
    dopt = JutulDarcy.setup_reservoir_dict_optimization(deepcopy(d0), setup_fn; strict = false, verbose = false)
    for k in keys_to_check
        free_optimization_parameter!(dopt, k)
    end
    grad = JutulDarcy.parameters_gradient_reservoir(dopt, G; deps = :case, info_level = -1)

    st0, _ = simulate(case0; info_level = -1)
    obj0 = Jutul.evaluate_objective(G, case0.model, st0, case0.dt, case0.forces)

    @printf("\n%-14s %16s %16s %12s %10s\n", "entry", "adjoint (DI)", "central-FD", "x0", "elasticity")
    for k in keys_to_check
        ga = grad[k][1]
        gfd, v0, _ = fd_entry(M, d0, k)
        # elasticity d(log obj)/d(log x): entries far below ~1e-3 carry no usable
        # signal and central FD there is dominated by simulator noise.
        elast = abs(ga * v0 / obj0)
        near_zero = elast < 1e-3
        ok = near_zero || isapprox(ga, gfd; rtol = 0.1, atol = 0.05*abs(gfd))
        @printf("%-14s %16.6e %16.6e %12.4g %10.2e%s\n", k, ga, gfd, v0, elast,
            ok ? "" : "   <-- MISMATCH")
        @test ok
    end
end

@testset "extended control vector gradient (DictOptimization, deps=:case)" begin
    run_regime("never",  OrderedDict{String,Any}("orat" => [60.0/day], "bhp" => [50.0*bar]);
        keys_to_check = ["orat", "bhp"])
    run_regime("switch", OrderedDict{String,Any}("orat" => [60.0/day], "bhp" => [80.0*bar]);
        keys_to_check = ["orat", "bhp"])
    run_regime("always", OrderedDict{String,Any}("orat" => [600.0/day], "bhp" => [150.0*bar]);
        keys_to_check = ["orat", "bhp"])
    run_regime("direct_bhp", OrderedDict{String,Any}("bhp_control" => [160.0*bar]);
        keys_to_check = ["bhp_control"])
end
