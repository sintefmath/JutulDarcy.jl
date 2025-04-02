@kwdef mutable struct ContractionFactorCuttingCriterion

    distance_function = r -> compute_distance(r)
    # Arrays for storing residuals/contraction factors
    history = nothing
    # Target number of nonlinear iterations for the timestep
    target_iterations = 5
    # Contraction factor parameters
    slow = 0.99
    fast = 0.1
    residual_history_length = 1
    # Violation counter and limits for action and timestep cut
    num_violations::Int     = 0
    num_violations_cut::Int = 3

end

function set_contraction_factor_cutting_criterion!(config; max_nonlinear_iterations = 5, kwargs...)

    cc = ContractionFactorCuttingCriterion(; kwargs...)
    config[:cutting_criterion] = cc
    config[:max_nonlinear_iterations] = max_nonlinear_iterations

end

function Jutul.cutting_criterion(cc::ContractionFactorCuttingCriterion, sim, dt, forces, it, max_iter, cfg, e, step_reports, relaxation)
    
    early_cut = false
    N = max(max_it - it + 1, 2)
    
    report = step_reports[end][:errors]
    dist = cc.distance_function(report)

    (it == 1) ? reset!(cc, dist, max_iter) : nothing

    # Store distance
    cc.history[:distance][it] = dist

    it0 = max(it - cc.residual_history_length, 1)
    Θ, Θ_target = compute_contraction_factor(cc.distances[it0:it], N)

    is_oscillating = oscillation(cc.history[:contraction_factor], cc.slow)
    
    good = all(Θ .<= max(Θ_target, cc.fast)) && !is_oscillating
    ok = all(Θ .<= cc.slow)
    bad = any(Θ .> cc.slow)

    if good
        # Convergence rate good, decrease number of violations
        cc.num_violations -= 1
        status = :good
    elseif ok
        # Convergence rate ok, keep number of violations
        status = :ok
    elseif bad
        # Not converging, increase number of violations
        cc.num_violations += 1
        status = :bad
    else
        # Sanity check - this should not happen
        error()
    end
    cc.history[:status][it] = status

    early_cut = cc.num_violations > cc.num_violations_cut
    step_reports[end][:cutting_criterion] = cc

    if cfg[:info_level] >= 2
        print_progress(cc, it)
    end

    return (relaxation, early_cut)

end

function reset!(cc::ContractionFactorCuttingCriterion, template, max_iter)

    cc.num_violations = 0

    history = Dict()
    history[:distance] = Array{typeof(template[1])}(undef, max_iter, length(template))
    history[:contraction_factor] = Array{typeof(template[1])}(undef, max_iter, length(template))
    history[:contraction_factor_target] = Array{typeof(template[1])}(undef, max_iter, length(template))
    history[:status] = Array{Symbol}(undef, max_iter, length(template))

    history[:contraction_factor][1] = NaN
    history[:contraction_factor_target][1] = NaN
    history[:status][1] = :none
    cc.history = history

end

function print_progress(cc::ContractionFactorCuttingCriterion, it, max_iter)

    round_local = x -> round(x; digits = 3)

    θ = cc.history[:contraction_factor][it]
    θ_target = cc.history[:contraction_factor_target][it]
    θ_slow = cc.slow
    θ_fast = cc.fast
    status = cc.history[:status][it]

    θ = round_local(θ)
    θ_target = round_local(θ_target)
    θ_slow = round_local(θ_slow)
    θ_fast = round_local(θ_fast)
    
    if status == :good
        θ_upper = max(θ_target, θ_fast)
    elseif status == :ok

    elseif status == :bad

    else
        error("Unknown status")
    end

    msg = "Convergence monitor (it = $it): "
    msg *= "status = $(cc.history[:status][it]), "
    msg *= "violations = $(cc.num_violations), "
    msg *= "violations = $(cc.num_violations), "
    msg *= "distance = $(cc.history[:distance][it]), "
    msg *= "contraction factor = $(cc.history[:contraction_factor][it])"
    Jutul.jutul_message("Convergence monitor", msg)

end