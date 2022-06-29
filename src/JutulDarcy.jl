module JutulDarcy
    import Jutul: number_of_cells, number_of_faces,
                  degrees_of_freedom_per_entity,
                  values_per_entity,
                  absolute_increment_limit, relative_increment_limit, maximum_value, minimum_value,
                  select_primary_variables_system!,
                  select_primary_variables_domain!,
                  initialize_primary_variable_ad!,
                  update_primary_variable!,
                  select_secondary_variables_flow_type!,
                  select_secondary_variables_system!,
                  select_secondary_variables_domain!,
                  select_secondary_variables_formulation!,
                  update_secondary_variable!,
                  default_value,
                  initialize_variable_value!,
                  initialize_variable_ad,
                  update_half_face_flux!,
                  update_accumulation!,
                  update_equation!,
                  minimum_output_variables,
                  select_equations,
                  select_equations_domain!,
                  select_equations_formulation!,
                  select_equations_system!,
                  select_output_variables,
                  setup_parameters,
                  count_entities,
                  count_active_entities,
                  associated_entity,
                  active_entities,
                  number_of_entities,
                  declare_entities,
                  get_neighborship

    import Jutul: setup_parameters_domain!, setup_parameters_system!, setup_parameters_context!, setup_parameters_formulation!
    import Jutul: fill_equation_entries!, update_linearized_system_equation!, check_convergence, update!, linear_operator, transfer, operator_nrows, matrix_layout, apply!
    import Jutul: apply_forces_to_equation!, convergence_criterion
    import Jutul: get_dependencies
    import Jutul: setup_forces, setup_state, setup_state!
    import Jutul: declare_pattern
    import Jutul: number_of_equations_per_entity
    using Jutul
    using ForwardDiff, StaticArrays, SparseArrays, LinearAlgebra, Statistics
    using AlgebraicMultigrid, Krylov
    # PVT
    using MultiComponentFlash
    using MAT
    using Tullio, LoopVectorization, Polyester, CUDA
    using TimerOutputs
    using DataStructures

    export reservoir_linsolve, get_1d_reservoir
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
    include("linsolve.jl")
    include("cpr.jl")
    include("deck_support.jl")
end # module
