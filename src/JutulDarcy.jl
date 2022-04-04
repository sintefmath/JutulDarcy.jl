module JutulDarcy
    import Jutul: number_of_cells, 
                  degrees_of_freedom_per_entity,
                  select_primary_variables_system!,
                  select_primary_variables_domain!,
                  select_secondary_variables_flow_type!,
                  select_secondary_variables_system!,
                  select_secondary_variables_domain!,
                  select_secondary_variables_formulation!,
                  update_secondary_variable!,
                  update_half_face_flux!,
                  update_accumulation!,
                  update_equation!,
                  minimum_output_variables,
                  select_equations,
                  select_equations_domain!,
                  select_equations_formulation!,
                  select_equations_system!,
                  select_output_variables,
                  count_entities,
                  count_active_entities,
                  declare_entities,
                  get_neighborship

    import Jutul: fill_equation_entries!, update_linearized_system_equation!, check_convergence
    using Jutul
    using ForwardDiff, StaticArrays, LinearAlgebra
    # PVT
    using MultiComponentFlash
    using MAT
    using Tullio, LoopVectorization, Polyester, CUDA
    using TimerOutputs
    include("types.jl")
    include("deck_types.jl")
    include("porousmedia_grids.jl")
    include("utils.jl")
    include("interpolation.jl")
    include("flux.jl")
    # Definitions for multiphase flow
    include("multiphase.jl")
    include("multiphase_secondary_variables.jl")
    # Compositional flow
    include("multicomponent/multicomponent.jl")

    # Blackoil
    include("blackoil/blackoil.jl")

    # Wells etc.
    include("facility/types.jl")

    include("facility/flux.jl")
    include("facility/wells.jl")
    include("facility/facility.jl")
    include("porousmedia.jl")
    # MRST inputs and test cases that use MRST input
    include("mrst_input.jl")
    # Various input tricks
    include("io.jl")
    include("cpr.jl")
    include("deck_support.jl")
end # module
