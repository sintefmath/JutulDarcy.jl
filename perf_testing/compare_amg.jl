# Compare the existing pressure-hierarchy update policies with the same case.
include("compare_solvers.jl")

policies = (
    "default" => NamedTuple(),
    "sparsity" => (update_type = :sparsity, amg_rebuild_tolerance = 0.0),
    "operators" => (update_type = :operators, amg_rebuild_tolerance = 0.0),
    "adaptive_20pct" => (amg_rebuild_tolerance = 0.2,),
    "no_partial" => (partial_update = false,),
)
selected = split(get(ENV, "PERF_AMG_POLICIES", join(first.(policies), " ")))
compare_solvers(; solver_policies = filter(p -> first(p) in selected, policies))
