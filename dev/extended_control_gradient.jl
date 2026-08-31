# Diagnostic: gradient w.r.t. an "extended control vector"
# =======================================================
#
# A producer is given a PRIMARY target (surface oil rate) AND a bhp lower
# limit. Both appear as separate entries in the flattened force vector:
# `target_Producer` and `limit_Producerbhp`.
#
# We compare the force adjoint (`solve_adjoint_forces`, `:all` facility
# variant) against central finite differences in four regimes:
#   never      - bhp limit never active
#   switch     - well drops from orat onto bhp mid-horizon
#   always     - well pinned to bhp limit the whole horizon
#   direct_bhp - well is BHP-controlled directly from forces (control baseline)
#
# Run:  julia --project=. dev/extended_control_gradient.jl
#
# FINDING (2026-08-31, Jutul 9riDr / JutulDarcy 0.3.10):
#   * `direct_bhp`: adjoint matches FD to ~6 sig figs -> the force adjoint is
#     correct for a DIRECT control target.
#   * `never/switch/always`: `limit_Producerbhp` adjoint entry is *identically
#     0.0* even when the bhp limit is the sole active control for 20/20 steps
#     (FD says it's ~5.3e-11, the same physical d(obj)/d(bhp) the `direct_bhp`
#     case gets right). And `target_Producer` gets a spurious ~0.45-0.53
#     regardless of whether the well is actually on its oil-rate target.
#   * Root cause: the adjoint re-evaluates the residual with a fresh
#     WellGroupConfiguration and re-runs `apply_well_limits!` at the converged
#     state, where `cond.bhp == bhp_limit` exactly, so the strict `<` test in
#     `check_well_limit` fails and NO switch is detected. The adjoint therefore
#     differentiates the (requested) oil-rate control equation on every step and
#     never sees the limit value in any residual.
#   The forward solve DOES store the correct operating control in
#   `state[:Facility][:WellGroupConfiguration].operating_controls` (that's what
#   `wells[:Producer][:control]` reads); the adjoint just doesn't use it.

using Jutul, JutulDarcy, LinearAlgebra, Printf

const bar = 1e5
const day = si_unit(:day)

function build_case(; orat_target = missing, bhp_floor = missing, bhp_control = missing,
                      nstep = 20, dtdays = 15.0)
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

    if !ismissing(bhp_control)
        # Well genuinely BHP-controlled from the forces (no limit switching).
        forces0 = setup_reservoir_forces(model;
            control = Dict(:Producer => ProducerControl(BottomHolePressureTarget(bhp_control))),
            set_default_limits = false)
    else
        forces0 = setup_reservoir_forces(model;
            control = Dict(:Producer => ProducerControl(SurfaceOilRateTarget(-orat_target))),
            limits = Dict(:Producer => (bhp = bhp_floor,)),
            set_default_limits = false)
    end
    forces = [forces0 for _ in 1:nstep]
    return JutulCase(model, dt, forces; state0 = state0, parameters = parameters)
end

function make_objective(case)
    t_tot = sum(case.dt)
    function G(model, state, dt, step_info, forces)
        orat = JutulDarcy.compute_well_qoi(model, state, forces, :Producer, SurfaceOilRateTarget)
        return (dt*orat)/t_tot
    end
    return G
end

function force_name_index(model, forces_single)
    _, cfg = Jutul.vectorize_forces(forces_single, model)
    idx = Dict{Symbol, Int}()
    gbase = 0
    for k in Jutul.submodels_symbols(model)
        haskey(cfg, k) || continue
        sub = cfg[k]
        lbase = 0
        i = 1
        for (_, m) in pairs(sub.meta)
            names = get(m, :names, Symbol[])
            for (j, nm) in enumerate(names)
                idx[nm] = gbase + lbase + j
            end
            lbase += sub.lengths[i]
            i += 1
        end
        gbase += sum(sub.lengths)
    end
    return idx
end

function adjoint_gradient(case)
    (; model, state0, parameters, dt, forces) = case
    G = make_objective(case)
    states, reports = simulate(case; output_substates = true, info_level = -1)
    _, _, grad_adj = Jutul.solve_adjoint_forces(model, states, reports, G, forces;
        state0 = state0, parameters = parameters)
    obj0 = Jutul.evaluate_objective(G, model, states, dt, forces)
    return grad_adj[1], obj0
end

function fd_entry(case, gidx; rel = 1e-3, absmin = 1e-8)
    (; model, state0, parameters, dt, forces) = case
    G = make_objective(case)
    f0 = forces[1]
    x0, cfg = Jutul.vectorize_forces(f0, model)
    eps = max(absmin, rel*abs(x0[gidx]))
    function ev(s)
        xp = copy(x0); xp[gidx] += s
        fp = Jutul.devectorize_forces(f0, model, xp, cfg)
        fvec = [fp for _ in dt]
        st, _ = simulate(state0, model, dt; parameters = parameters, forces = fvec, info_level = -1)
        return Jutul.evaluate_objective(G, model, st, dt, fvec)
    end
    return (ev(eps) - ev(-eps))/(2*eps), x0[gidx], eps
end

regimes = (
    never      = (orat_target =  60.0/day, bhp_floor =  70*bar),
    switch     = (orat_target =  60.0/day, bhp_floor = 185*bar),
    always     = (orat_target = 600.0/day, bhp_floor = 150*bar),
    direct_bhp = (bhp_control = 160*bar,),
)

for (rname, p) in pairs(regimes)
    println("\n" * "="^70)
    println("REGIME $rname   $(p)")
    println("="^70)
    case = build_case(; p...)

    res = simulate_reservoir(case; info_level = -1)
    ctrl = res.wells[:Producer][:control]
    orat = res.wells[:Producer][:orat]
    bhp  = res.wells[:Producer][:bhp]
    @printf("control per step : %s\n", string(ctrl))
    @printf("orat  (m3/day)   : %s\n", string(round.(orat .* day, digits = 1)))
    @printf("bhp   (bar)      : %s\n", string(round.(bhp ./ bar, digits = 1)))

    idx = force_name_index(case.model, case.forces[1])
    g, obj = adjoint_gradient(case)
    @printf("\nobjective = %.6e   (mean produced oil rate, m3/s)\n", obj)
    @printf("%-22s %14s %14s %10s\n", "entry", "adjoint", "central-FD", "x0")
    for nm in sort(collect(keys(idx)))
        gi = idx[nm]
        ga = g[gi]
        gfd, x0v, ep = fd_entry(case, gi)
        flag = ""
        denom = max(abs(ga), abs(gfd), 1e-300)
        if abs(ga - gfd) > 0.05*denom + 1e-12
            flag = "  <-- MISMATCH"
        end
        @printf("%-22s %14.6e %14.6e %10.3g%s\n", nm, ga, gfd, x0v, flag)
    end
end
