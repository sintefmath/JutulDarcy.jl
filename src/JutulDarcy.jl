__precompile__(false)

"""
$(README)

---
##  Module exports:

$(EXPORTS)

"""
module JutulDarcy
    import Jutul:
        number_of_cells, number_of_faces,
        degrees_of_freedom_per_entity,
        values_per_entity,
        absolute_increment_limit, relative_increment_limit, maximum_value, minimum_value,
        select_primary_variables!,
        select_secondary_variables!,
        select_equations!,
        select_parameters!,
        select_minimum_output_variables!,
        initialize_primary_variable_ad!,
        update_primary_variable!,
        update_secondary_variable!,
        default_value,
        initialize_variable_value!,
        initialize_variable_value,
        initialize_variable_ad!,
        update_half_face_flux!,
        update_accumulation!,
        update_equation!,
        setup_parameters,
        count_entities,
        count_active_entities,
        associated_entity,
        active_entities,
        number_of_entities,
        declare_entities,
        get_neighborship

    import Jutul: update_preconditioner!, partial_update_preconditioner!

    import Jutul: fill_equation_entries!, update_linearized_system_equation!, check_convergence, update!, linear_operator, transfer, operator_nrows, matrix_layout, apply!
    import Jutul: apply_forces_to_equation!, convergence_criterion
    import Jutul: get_dependencies
    import Jutul: setup_forces, setup_state, setup_state!
    import Jutul: declare_pattern
    import Jutul: number_of_equations_per_entity
    import Jutul: update_equation_in_entity!, update_cross_term_in_entity!, local_discretization, discretization
    import Jutul: FiniteVolumeGlobalMap, TrivialGlobalMap
    import Jutul: @tic

    using Jutul
    using ForwardDiff, StaticArrays, SparseArrays, LinearAlgebra, Statistics
    using AlgebraicMultigrid
    # PVT
    using MultiComponentFlash
    using MAT
    using Tullio, LoopVectorization, Polyester
    using TimerOutputs
    using PrecompileTools
    using Dates
    import DataStructures: OrderedDict
    using DocStringExtensions


    include("types.jl")
    include("deck_types.jl")
    include("porousmedia_grids.jl")
    include("utils.jl")
    include("interpolation.jl")
    include("flux.jl")
    # Definitions for multiphase flow
    include("multiphase.jl")
    include("variables/variables.jl")
    # Compositional flow
    include("multicomponent/multicomponent.jl")

    # Blackoil
    include("blackoil/blackoil.jl")

    include("thermal/thermal.jl")

    # Wells etc.
    include("facility/facility.jl")

    include("porousmedia.jl")
    # MRST inputs and test cases that use MRST input
    # and .DATA file simulation
    include("input_simulation/input_simulation.jl")
    # Initialization by equilibriation
    include("init/init.jl")
    # Corner point grids
    include("cpgrid/cpgrid.jl")
    # Gradients, objective functions, etc
    include("gradients/gradients.jl")

    # Various input tricks
    include("io.jl")
    include("linsolve.jl")
    include("cpr.jl")
    include("deck_support.jl")
    include("regions/regions.jl")
    include("test_utils/test_utils.jl")
    include("forces/forces.jl")

    include("formulations/formulations.jl")

    include("ext.jl")

    include("InputParser/InputParser.jl")
    using .InputParser
    import .InputParser: parse_data_file

    export
        reservoir_linsolve,
        parse_data_file,
        get_1d_reservoir,
        DeckViscosity,
        DeckShrinkage,
        CPRPreconditioner,
        MuBTable,
        ConstMuBTable,
        DeckDensity, RelativePermeabilities,
        ThreePhaseCompositionalDensitiesLV,
        PhaseMassFractions,
        PhaseMassFractions,
        ThreePhaseLBCViscositiesLV,
        plot_well!,
        plot_well_results,
        plot_reservoir_simulation_result,
        plot_reservoir,
        simulate_reservoir_parray,
        setup_reservoir_simulator_parray,
        component_mass_fluxes!, update_total_masses!,
        table_to_relperm,
        MultiPhaseSystem, ImmiscibleSystem, SinglePhaseSystem,
        AqueousPhase, LiquidPhase, VaporPhase,
        number_of_phases, get_short_name, get_name,
        SourceTerm,
        Pressure, Saturations, TotalMasses, TotalMass,
        fluid_volume, pore_volume,
        MinimalTPFAGrid,
        compute_peaceman_index,
        discretized_domain_tpfv_flow,
        discretized_domain_well,
        StandardBlackOilSystem,
        PhaseRelPerm,
        FlowBoundaryCondition,
        ReservoirSimResult,
        reservoir_domain,
        reservoir_model,
        setup_reservoir_model,
        setup_reservoir_simulator,
        simulate_reservoir,
        setup_reservoir_state,
        setup_reservoir_forces,
        full_well_outputs,
        well_output,
        well_symbols,
        wellgroup_symbols,
        available_well_targets,
        BlackOilUnknown,
        BlackOilX,
        ReservoirFromWellFlowCT,
        FacilityFromWellFlowCT,
        WellFromFacilityFlowCT,
        TotalSurfaceMassRate,
        WellGroup,
        DisabledControl,
        HistoryMode, PredictionMode,
        Wells,
        TotalMassVelocityMassFractionsFlow
        BottomHolePressureTarget, TotalRateTarget, SinglePhaseRateTarget, DisabledTarget
        SurfaceLiquidRateTarget, SurfaceOilRateTarget, SurfaceWaterRateTarget, SurfaceGasRateTarget,
        PerforationMask,
        WellDomain, MultiSegmentWell,
        TotalMassFlux, PotentialDropBalanceWell, SegmentWellBoreFrictionHB,
        InjectorControl, ProducerControl, SinglePhaseRateTarget, BottomHolePressureTarget,
        Perforations,
        MixedWellSegmentFlow,
        segment_pressure_drop,
        setup_well, setup_vertical_well,
        well_mismatch,
        simulate_data_file, setup_case_from_data_file,
        get_test_setup, get_well_from_mrst_data,
        setup_case_from_mrst,
        simulate_mrst_case,
        MultiPhaseCompositionalSystemLV,
        StandardVolumeSource, VolumeSource, MassSource,
        OverallMoleFractions,
        ImmiscibleSaturation,
        ThermalSystem,
        PhaseMassDensities, ConstantCompressibilityDensities,
        BrooksCoreyRelPerm, TabulatedRelPermSimple

    @compile_workload begin
        precompile_darcy_multimodels()
        # We run a tiny MRST case to precompile the .MAT file loading
        spe1_path = joinpath(pathof(JutulDarcy), "..", "..", "test", "mrst", "spe1.mat")
        if isfile(spe1_path)
            simulate_mrst_case(spe1_path, info_level = -1, verbose = false)
        end
    end
end # module
