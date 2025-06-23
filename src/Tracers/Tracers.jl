module Tracers
    using Jutul, JutulDarcy
    import JutulDarcy: AbstractPhase, get_phases, darcy_phase_mass_flux, phase_upwind, MultiPhaseSystem
    import StaticArrays: setindex

    export SinglePhaseTracer, MultiPhaseTracer, add_tracers_to_model!, number_of_tracers

    include("tracer_variants.jl")
    include("variables.jl")
    include("equations.jl")
    include("facility.jl")
    include("api.jl")
    include("polymer.jl")
end
