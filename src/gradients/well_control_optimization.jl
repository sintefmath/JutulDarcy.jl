# Well control optimization on top of Jutul's DictOptimization
# ============================================================
#
# Workflow (MRST `optimizeEggEnsemble` style, minus the ensemble):
#
#   1. Start from a `JutulCase` with a base schedule.
#   2. Declare control periods explicitly as groups of report steps.
#   3. Declare, per well, exactly which quantities are free (primary target
#      and/or limits) and their box bounds (absolute and/or relative).
#   4. `optimize_well_controls` runs the adjoint-based quasi-Newton optimizer
#      (`Jutul.unit_box_bfgs` via the `:lbfgs` backend).
#
# A well quantity is treated as the *target* when it matches the well's base
# control target, otherwise as a *limit* (it is written into
# `forces[:Facility].limits[well]`). The adjoint follows whichever control is
# active during the simulation (see the reference-mode handling in
# `src/facility/controls.jl`).

"""
Symbol -> well target constructor for the quantities that can be optimized.
"""
const WELL_CONTROL_TARGET_TYPES = Dict{Symbol, Any}(
    :bhp  => BottomHolePressureTarget,
    :rate => TotalRateTarget,
    :orat => SurfaceOilRateTarget,
    :wrat => SurfaceWaterRateTarget,
    :grat => SurfaceGasRateTarget,
    :lrat => SurfaceLiquidRateTarget,
)

struct WellControlDOF
    well::Symbol
    quantity::Symbol
    is_target::Bool          # true: primary control target, false: operating limit
    periods::Vector{Int}     # which control-period indices this DOF is active in
    constant::Bool           # one shared value across all its periods
    initial::Vector{Float64} # magnitude, per active period (already expanded)
    abs_min::Float64
    abs_max::Float64
    rel_min::Float64
    rel_max::Float64
    scaler::Union{Symbol, Missing}
end

struct WellControlOptimization
    dopt::Jutul.DictParameters
    base_case::JutulCase
    control_periods::Vector{Vector{Int}}
    step_to_period::Vector{Int}   # length nstep, 0 == not under optimization control
    dofs::Vector{WellControlDOF}
    objective
    setup_function
end

(copt::WellControlOptimization)(prm::AbstractDict) = copt.setup_function(prm, missing)

function Base.show(io::IO, ::MIME"text/plain", copt::WellControlOptimization)
    np = length(copt.control_periods)
    println(io, "WellControlOptimization: $(length(copt.dofs)) DOF specs over $np control period(s)")
    for (i, p) in enumerate(copt.control_periods)
        println(io, "  period $i: steps $(_compact_range(p))")
    end
    for d in copt.dofs
        kind = d.is_target ? "target" : "limit"
        rng = d.constant ? "constant" : "periods $(d.periods)"
        println(io, "  $(d.well).$(d.quantity) ($kind, $rng)  box [$(d.abs_min), $(d.abs_max)]" *
            (d.scaler === missing ? "" : "  scaler=$(d.scaler)"))
    end
end

function _compact_range(v)
    isempty(v) && return "[]"
    s = sort(collect(v))
    (s == collect(first(s):last(s))) && return "$(first(s)):$(last(s))"
    return string(s)
end

"""
    control_is_injector(base_forces, well) -> Bool

Look up whether `well` is an injector in the base schedule.
"""
function control_is_injector(base_forces, well)
    ctrl = base_forces[:Facility].control[well]
    return ctrl isa InjectorControl
end

"""
Signed physical value from an optimization magnitude, following JutulDarcy's
sign conventions (injection rates positive, production rates negative, pressures
positive).
"""
function wc_signed_value(quantity::Symbol, magnitude, is_injector::Bool)
    if quantity == :bhp
        return magnitude
    else
        return is_injector ? magnitude : -magnitude
    end
end

function _wc_base_forces_vector(case::JutulCase)
    f = case.forces
    n = length(case.dt)
    if f isa AbstractVector
        length(f) == n || error("case.forces has $(length(f)) entries but there are $n report steps.")
        return [deepcopy(f[i]) for i in 1:n]
    else
        return [deepcopy(f) for _ in 1:n]
    end
end

function _wc_normalize_periods(control_periods, nstep)
    periods = Vector{Vector{Int}}()
    for p in control_periods
        pv = sort(collect(Int, p))
        isempty(pv) && error("Empty control period.")
        (first(pv) >= 1 && last(pv) <= nstep) ||
            error("Control period $pv out of range 1:$nstep.")
        push!(periods, pv)
    end
    step_to_period = zeros(Int, nstep)
    for (i, pv) in enumerate(periods)
        for s in pv
            step_to_period[s] == 0 ||
                error("Report step $s appears in more than one control period.")
            step_to_period[s] = i
        end
    end
    return periods, step_to_period
end

function _wc_expand_initial(spec_initial, base_magnitude, active_periods)
    npa = length(active_periods)
    if spec_initial === missing
        base_magnitude === missing &&
            error("No initial value given and none found in the base schedule.")
        return fill(Float64(base_magnitude), npa)
    elseif spec_initial isa Number
        return fill(Float64(spec_initial), npa)
    else
        v = collect(Float64, spec_initial)
        length(v) == npa || error("initial has length $(length(v)), expected $npa (active periods).")
        return v
    end
end

function _wc_base_magnitude(base_forces_vec, periods, well, quantity, is_target)
    # value from the first step of each active period; use the first period's
    # value as the representative base magnitude for relative bounds
    p1 = first(periods)
    bf = base_forces_vec[first(p1)]
    fac = bf[:Facility]
    if is_target
        return abs(fac.control[well].target.value)
    else
        lim = get(fac.limits, well, nothing)
        if lim === nothing || !haskey(lim, quantity)
            return missing
        end
        return abs(lim[quantity])
    end
end

"""
    setup_well_control_optimization(case, control_periods, controls; kwargs...)

Set up a well-control optimization problem.

# Arguments
- `case::JutulCase`: base case (model, base schedule, `state0`, parameters).
- `control_periods`: vector of report-step groups, e.g. `[1:20, 21:40]`. Groups
  must be disjoint. Report steps not covered by any group keep their base forces
  and are not differentiated with respect to controls.
- `controls`: vector of `NamedTuple` DOF specs, one per (well, quantity):
    - `well::Symbol`, `quantity::Symbol` (`:bhp`, `:rate`, `:orat`, `:wrat`,
      `:grat`, `:lrat`) — required.
    - `abs_min`, `abs_max`: absolute bounds on the magnitude (physical units).
    - `rel_min`, `rel_max`: bounds relative to the base magnitude.
    - `initial`: starting magnitude (scalar or one per active period). Defaults
      to the base schedule value; required for a limit absent from the base
      schedule.
    - `constant::Bool = false`: use a single shared value across the periods.
    - `periods`: indices into `control_periods` the DOF is active in (default all).
    - `scaler::Symbol = :linear_limits`: DictOptimization scaler.

# Keyword arguments
- `objective = :npv`: `:npv` (uses [`npv_objective`](@ref)) or a function
  `(model, state, dt, step_info, forces) -> Float64`.
- `npv_arg = NamedTuple()`: forwarded to `npv_objective` when `objective = :npv`.
- `strict`, `verbose`: forwarded to `DictParameters`.

Returns a [`WellControlOptimization`](@ref). Call
[`optimize_well_controls`](@ref) to run the optimizer, then call the returned
object on the optimized parameter dict to rebuild the `JutulCase`.
"""
function setup_well_control_optimization(case::JutulCase, control_periods, controls;
        objective = :npv,
        npv_arg = NamedTuple(),
        strict = true,
        verbose = true,
    )
    model = case.model
    nstep = length(case.dt)
    base_forces_vec = _wc_base_forces_vector(case)
    periods, step_to_period = _wc_normalize_periods(control_periods, nstep)
    np = length(periods)

    wells = Set(well_symbols(model))
    dofs = WellControlDOF[]
    for (i, c) in enumerate(controls)
        haskey(c, :well) && haskey(c, :quantity) ||
            error("controls[$i] must have :well and :quantity.")
        well = c.well::Symbol
        quantity = c.quantity::Symbol
        well in wells || error("controls[$i]: unknown well $well.")
        haskey(WELL_CONTROL_TARGET_TYPES, quantity) ||
            error("controls[$i]: unsupported quantity $quantity " *
                  "(supported: $(sort(collect(keys(WELL_CONTROL_TARGET_TYPES))))).")

        active_periods = haskey(c, :periods) ? sort(collect(Int, c.periods)) : collect(1:np)
        all(1 .<= active_periods .<= np) ||
            error("controls[$i]: periods $(active_periods) out of range 1:$np.")

        base_ctrl = base_forces_vec[first(first(periods))][:Facility].control[well]
        base_ctrl isa DisabledControl &&
            error("controls[$i]: well $well is disabled in the base schedule.")
        primary_sym = translate_target_to_symbol(base_ctrl.target)
        is_target = quantity == primary_sym

        base_mag = _wc_base_magnitude(base_forces_vec, periods, well, quantity, is_target)
        initial = _wc_expand_initial(get(c, :initial, missing), base_mag, active_periods)

        abs_min = Float64(get(c, :abs_min, 0.0))
        abs_max = Float64(get(c, :abs_max, Inf))
        rel_min = Float64(get(c, :rel_min, -Inf))
        rel_max = Float64(get(c, :rel_max, Inf))
        scaler = get(c, :scaler, :linear_limits)
        constant = Bool(get(c, :constant, false))

        # Resolve relative bounds against the base magnitude (or the initial if
        # there is no base magnitude).
        ref = base_mag === missing ? first(initial) : base_mag
        lo = max(abs_min, isfinite(rel_min) ? rel_min*ref : -Inf)
        hi = min(abs_max, isfinite(rel_max) ? rel_max*ref : Inf)
        isfinite(hi) || error("controls[$i]: no finite upper bound (set abs_max or rel_max).")
        lo < hi || error("controls[$i]: empty box [$lo, $hi].")
        all(lo .<= initial .<= hi) ||
            @warn "controls[$i] ($well.$quantity): initial $(initial) outside box [$lo, $hi]; will be clamped."

        push!(dofs, WellControlDOF(well, quantity, is_target, active_periods, constant,
            clamp.(initial, lo, hi), lo, hi, rel_min, rel_max, scaler))
    end

    # Build the optimization dict: well => quantity => per-period magnitudes.
    dict = OrderedDict{String, Any}()
    for d in dofs
        wkey = String(d.well)
        haskey(dict, wkey) || (dict[wkey] = OrderedDict{String, Any}())
        dict[wkey][String(d.quantity)] = d.constant ? [first(d.initial)] : copy(d.initial)
    end

    F = (prm, step_info = missing) -> wc_rebuild_case(prm, case, base_forces_vec,
        dofs, periods, step_to_period)

    dopt = setup_reservoir_dict_optimization(dict, F; strict = strict, verbose = verbose)
    for d in dofs
        name = [String(d.well), String(d.quantity)]
        free_optimization_parameter!(dopt, name;
            abs_min = d.abs_min,
            abs_max = d.abs_max,
            scaler = d.scaler,
        )
    end

    obj = objective === :npv ? wc_npv_objective(case, base_forces_vec; npv_arg...) : objective

    return WellControlOptimization(dopt, case, periods, step_to_period, dofs, obj, F)
end

"""
Rebuild a `JutulCase` from an optimization parameter dict, overriding only the
controlled wells on the controlled steps. All steps within one control period
share a single forces object so the adjoint can reuse sparsity per unique force.
"""
function wc_rebuild_case(prm::AbstractDict, case::JutulCase, base_forces_vec,
        dofs::Vector{WellControlDOF}, periods, step_to_period)
    model = case.model
    nstep = length(case.dt)
    # One representative force object per period (deep-copied from the base).
    period_force = Dict{Int, Any}()
    for (pi, pv) in enumerate(periods)
        period_force[pi] = deepcopy(base_forces_vec[first(pv)])
    end
    for d in dofs
        wkey = String(d.well)
        qkey = String(d.quantity)
        vals = prm[wkey][qkey]
        is_inj = control_is_injector(base_forces_vec[first(first(periods))], d.well)
        TT = WELL_CONTROL_TARGET_TYPES[d.quantity]
        for (k, pi) in enumerate(d.periods)
            mag = d.constant ? vals[1] : vals[k]
            sval = wc_signed_value(d.quantity, mag, is_inj)
            f = period_force[pi]
            fac = f[:Facility]
            if d.is_target
                ctrl = fac.control[d.well]
                new_ctrl = replace_target(ctrl, TT(sval))
                fac.control[d.well] = new_ctrl
                lim = get(fac.limits, d.well, nothing)
                fac.limits[d.well] = lim === nothing ? as_limit(new_ctrl.target) :
                    merge(lim, as_limit(new_ctrl.target))
            else
                lim = get(fac.limits, d.well, nothing)
                add = NamedTuple{(d.quantity,)}((sval,))
                fac.limits[d.well] = lim === nothing ? add : merge(lim, add)
            end
        end
    end
    forces = Vector{Any}(undef, nstep)
    for s in 1:nstep
        p = step_to_period[s]
        forces[s] = p == 0 ? base_forces_vec[s] : period_force[p]
    end
    forces = identity.(forces)
    return JutulCase(model, case.dt, forces; state0 = case.state0, parameters = case.parameters)
end

"""
    wc_npv_objective(case, base_forces_vec; kwarg...)

Build an NPV objective closure for [`setup_well_control_optimization`](@ref),
picking up injector / producer lists from the base schedule. Keyword arguments
are forwarded to [`npv_objective`](@ref) (`oil_price`, `discount_rate`, ...).
"""
function wc_npv_objective(case::JutulCase, base_forces_vec; kwarg...)
    ctrls = base_forces_vec[1][:Facility].control
    injectors = Symbol[w for (w, c) in pairs(ctrls) if c isa InjectorControl]
    producers = Symbol[w for (w, c) in pairs(ctrls) if c isa ProducerControl]
    dt = case.dt
    return (model, state, t, step_info, forces) -> npv_objective(model, state, t, step_info, forces;
        injectors = injectors,
        producers = producers,
        timesteps = dt,
        kwarg...)
end

"""
    well_control_optimization_problem(copt::WellControlOptimization; deps = :case, kwarg...)

Standalone `Jutul.DictOptimization.JutulOptimizationProblem` for a well-control
problem: `f, g = prob(x)` evaluates the NPV-style objective and its adjoint
gradient at a (scaled) control vector `x`, `prob.x0` is the scaled initial
guess, and `prob.descale` maps back to physical magnitudes. Useful for gradient
checks and plugging into an external optimizer.
"""
# Dense differentiation over the (small) control vector. The sparse pattern is
# detected at the FIRST gradient evaluation and then frozen: a well limit that
# is inactive there would keep a hard-zero gradient for the rest of the
# optimization even after it becomes active (verified in
# dev/kink_cache_isolation.jl / dev/kink_ministep_anatomy.jl). Dense mode does
# not have this failure and costs little for typical control-DOF counts.
const WELL_CONTROL_BACKEND_ARG = (
    use_sparsity = false,
    di_sparse = false,
    single_step_sparsity = false,
    do_prep = true,
)

function well_control_optimization_problem(copt::WellControlOptimization;
        deps = :case,
        simulator_arg = (info_level = -1, end_report = false),
        backend_arg = WELL_CONTROL_BACKEND_ARG,
        kwarg...)
    # output_substates is always enforced so that every ministep enters the
    # adjoint solve (control/limit switching can happen within a report step).
    simulator_arg = merge(simulator_arg, (output_substates = true,))
    sim, cfg = setup_simulator_for_reservoir_optimization(copt.dopt, copt.setup_function,
        missing, missing, simulator_arg)
    return Jutul.DictOptimization.optimization_problem(copt.dopt, copt.objective, copt.setup_function;
        deps = deps, simulator = sim, config = cfg, backend_arg = backend_arg, kwarg...)
end

"""
    optimize_well_controls(copt::WellControlOptimization; maximize = true, kwarg...)

Run the optimizer on a problem from [`setup_well_control_optimization`](@ref).
Returns the optimized parameter dict; call `copt(prm)` to get the tuned
`JutulCase`. Keyword arguments are forwarded to
[`optimize_reservoir`](@ref) / `Jutul.optimize` (`max_it`, `grad_tol`,
`lin_eq`, ...).
"""
function optimize_well_controls(copt::WellControlOptimization;
        maximize = true,
        optimizer = :lbfgs,
        backend_arg = WELL_CONTROL_BACKEND_ARG,
        kwarg...)
    return optimize_reservoir(copt.dopt, copt.objective;
        deps = :case,
        maximize = maximize,
        optimizer = optimizer,
        backend_arg = backend_arg,
        kwarg...)
end
