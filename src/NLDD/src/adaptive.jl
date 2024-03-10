"""
    active_status(strategy::AbstractNLDDStrategy, log, iteration)

Determine active status (= use ASPEN/NLDD) for given iteration and log
"""
function active_status(strategy::AbstractNLDDStrategy, log, iteration)
    return true
end

"""
    active_status_first_iteration(strategy::AbstractNLDDStrategy, log)

Special version called for first iteration (start up)
"""
function active_status_first_iteration(strategy::AbstractNLDDStrategy, log)
    return true
end

"""
    active_status(strategy::DefaultNLDDStrategy, log, iteration)

Default strategy, just check `max_local_iterations` and disable advanced
algorithm if we have exceeded the limit for this step.
"""
function active_status(strategy::DefaultNLDDStrategy, log, iteration)
    # We should have exactly iteration-1 reports available (since we are solving
    # the current iteration)
    @assert length(log.reports) == iteration-1
    return iteration < strategy.max_local_iterations
end

function active_status_first_iteration(strategy::DistanceNLDDStrategy, log)
    start_active = !strategy.start_newton
    return start_active || log.forces_changed
end

function active_status(strategy::DistanceNLDDStrategy, log, it)
    inc_step = strategy.include_step
    inc_hist = strategy.include_history
    active = log.active
    completed_its = it - 1
    if completed_its > 1
        # We only know distances once two iterations have been
        # performed
        reports = log.reports
        d_first = distance_to_convergence(first(reports))
        d_prev = distance_to_convergence(reports[end-1])
        d = distance_to_convergence(reports[end])
        # Budget approach
        budget = strategy.iteration_budget
        budget_exceeded = false
        if isfinite(budget)
            budget_exceeded_hist = d > d_first*(1.0 - completed_its/budget)
            budget_exceeded_step = (d - d_prev) < d_first/budget
            budget_exceeded = (budget_exceeded_hist && inc_hist) || 
                              (budget_exceeded_step && inc_step)
        end
        # Rate approach
        rate = strategy.reduction_rate
        too_slow = false
        if isfinite(rate)
            too_slow_hist = d > d_first - rate*completed_its
            too_slow_step = (d - d_prev) < rate
            budget_exceeded = (too_slow_hist && inc_hist) || 
                              (too_slow_step && inc_step)
        end
        active = budget_exceeded || too_slow
        # @info "Switching $it" budget_exceeded too_slow
    else
        active = false
        # @info "Keeping $active $it"
    end
    # active = false
    return active
end

function distance_to_convergence(report)
    errors = report[:errors]
    if haskey(errors, :Reservoir)
        # Multi model
        # active = [:Reservoir]
        active = keys(errors)
        cnv = 0
        mb = 0
        for ek in active
            e = errors[ek]
            c = first(e).criterions
            if haskey(c, :CNV)
                cnv = max(maximum(c.CNV.errors), cnv)
                mb = max(maximum(c.MB.errors), mb)
            end
        end
    else
        # Scalar model
        c = first(errors.criterions)
        mb = maximum(c.MB.errors)
        cnv = maximum(c.CNV.errors)
    end
    return distance_to_convergence(mb, cnv)
end

function distance_to_convergence(mb_max, cnv_max)
    # TODO: These shouldn't be hard-coded.
    tol_cnv = 1e-3
    tol_mb = 1e-7

    M = log10(mb_max) - log10(tol_mb)
    C = log10(cnv_max) - log10(tol_cnv)

    return sqrt(max(M, 0)^2 + max(C, 0)^2)
end

function active_status_first_iteration(strategy::AlternatingNLDDStrategy, log)
    return !(strategy.start_newton)
end

function active_status(strategy::AlternatingNLDDStrategy, log, it)
    N = it + Int64(strategy.start_newton)
    return mod(N, 2) == 1
end

function active_status_first_iteration(strategy::OnlyWellsNLDDStrategy, log)
    return log.forces_changed
end

function active_status(strategy::OnlyWellsNLDDStrategy, log, it)
    return log.forces_changed
end


# Delta version
function active_status_first_iteration(strategy::VariableDeltaNLDDStrategy, log)
    return !(strategy.start_newton)
end

function active_status(strategy::VariableDeltaNLDDStrategy, log, iteration)
    # We should have exactly iteration-1 reports available (since we are solving
    # the current iteration)
    @assert length(log.reports) == iteration-1
    report = log.reports[end][:update]
    if haskey(report, :Reservoir)
        report = report[:Reservoir]
    end
    active = false
    for (k, v) in strategy.tolerances
        val, t = v
        if haskey(report, k)
            actual = report[k][t]
            active = active || actual > val
        end
    end
    return active
end

# Hybrid version
function active_status_first_iteration(strategy::HybridNLDDstrategy, log)
    return active_status_first_iteration(strategy.activation_strategy, log)
end

function active_status(strategy::HybridNLDDstrategy, log, iteration)
    if log.active
        active = active_status(strategy.disable_strategy, log, iteration)
    else
        active = active_status(strategy.activation_strategy, log, iteration)
    end
    return active
end
