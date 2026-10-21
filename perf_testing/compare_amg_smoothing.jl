# Reuse the full-simulation driver to test smoothing strength with existing APIs.
include("compare_solvers.jl")
policies = (
    "gs" => (amg_arg = (smoother_type = :gauss_seidel,),),
    "pressure_ilu" => (amg_arg = (smoother_type = :ilu0,),),
    "spai2" => (amg_arg = (smoother_type = :spai0, npre = 2),),
    "spai4" => (amg_arg = (smoother_type = :spai0, npre = 4),),
)
selected = split(get(ENV, "PERF_AMG_POLICIES", join(first.(policies), " ")))
compare_solvers(; solver_policies = filter(p -> first(p) in selected, policies))
