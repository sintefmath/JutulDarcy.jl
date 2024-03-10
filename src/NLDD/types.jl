struct NLDDSimulator <: JutulSimulator
    simulator
    executor
    partition
    subdomain_neighbors
    subdomain_simulators
    storage
end

mutable struct NLDDSolveLog
    active
    dt
    reports
    forces_changed
    forces
end

function NLDDSolveLog(active = true, dt = NaN; reports = [])
    NLDDSolveLog(active, dt, reports, false, nothing)
end

abstract type AbstractNLDDStrategy end

export DefaultNLDDStrategy, DistanceNLDDStrategy,
       AlternatingNLDDStrategy, OnlyWellsNLDDStrategy,
       VariableDeltaNLDDStrategy

struct DefaultNLDDStrategy <: AbstractNLDDStrategy
    max_local_iterations
    function DefaultNLDDStrategy(max_its = Inf)
        new(max_its)
    end
end

struct DistanceNLDDStrategy <: AbstractNLDDStrategy
    iteration_budget::Union{Int64, Float64} # Allow Inf
    reduction_rate::Float64
    start_newton::Bool
    include_step::Bool
    include_history::Bool
    function DistanceNLDDStrategy(N_opt = 6, red_rate = 1.0, newt = true; history = true, step = true)
        new(N_opt, red_rate, newt, step, history)
    end
end

struct AlternatingNLDDStrategy <: AbstractNLDDStrategy
    start_newton::Bool
    function AlternatingNLDDStrategy(newt = true)
        new(newt)
    end
end

struct OnlyWellsNLDDStrategy <: AbstractNLDDStrategy
end

struct VariableDeltaNLDDStrategy <: AbstractNLDDStrategy
    start_newton::Bool
    tolerances
    function VariableDeltaNLDDStrategy(newt = true; tol...)
        tol_vec = []
        for (k, v) in tol
            push!(tol_vec, (k, v))
        end
        new(newt, tol_vec)
    end
end

struct HybridNLDDstrategy <: AbstractNLDDStrategy
    activation_strategy
    disable_strategy
end
