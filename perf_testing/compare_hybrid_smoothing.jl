include("compare_solvers.jl")
const KP = Jutul.KAPreconditioners
Base.include(KP, joinpath(@__DIR__, "hybrid_smoother_candidate.jl"))
policies = (
    "gs" => (amg_arg = (smoother_type = :gauss_seidel,),),
    "spai1" => (amg_arg = (smoother_type = :spai0,),),
    "hybrid8k" => (amg_arg = (smoother = KP.PerfHybridSmoother(8192),),),
    "hybrid32k" => (amg_arg = (smoother = KP.PerfHybridSmoother(32768),),),
    "hybrid32k_2" => (amg_arg = (smoother = KP.PerfHybridSmoother(32768; steps = 2),),),
)
selected = split(get(ENV, "PERF_AMG_POLICIES", join(first.(policies), " ")))
compare_solvers(; solver_policies = filter(p -> first(p) in selected, policies))
