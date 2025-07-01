"""
$(README)
"""
module JutulDarcy
    export MultiPhaseSystem, ImmiscibleSystem, SinglePhaseSystem
    export reservoir_linsolve
    export get_1d_reservoir
    export DeckPhaseViscosities
    export DeckShrinkageFactors
    export CPRPreconditioner
    export MuBTable
    export ConstMuBTable
    export DeckPhaseMassDensities, RelativePermeabilities
    export ThreePhaseCompositionalDensitiesLV
    export PhaseMassFractions
    export PhaseMassFractions
    export ThreePhaseLBCViscositiesLV
    export plot_well!
    export plot_well_results
    export plot_reservoir_simulation_result
    export plot_reservoir
    export simulate_reservoir_parray
    export setup_reservoir_simulator_parray
    export component_mass_fluxes!, update_total_masses!
    export table_to_relperm
    export AqueousPhase, LiquidPhase, VaporPhase
    export number_of_phases
    export SourceTerm
    export Pressure, Saturations, TotalMasses, TotalMass
    export fluid_volume, pore_volume
    export compute_peaceman_index
    export discretized_domain_tpfv_flow
    export discretized_domain_well
    export StandardBlackOilSystem
    export PhaseRelativePermeability
    export FlowBoundaryCondition
    export ReservoirSimResult
    export reservoir_domain
    export reservoir_model
    export setup_reservoir_model
    export setup_reservoir_simulator
    export simulate_reservoir
    export setup_reservoir_state
    export setup_reservoir_forces
    export full_well_outputs
    export well_output
    export well_symbols
    export wellgroup_symbols
    export available_well_targets
    export BlackOilUnknown
    export BlackOilX
    export TotalSurfaceMassRate
    export WellGroup
    export DisabledControl
    export Wells
    export TotalMassVelocityMassFractionsFlow
    export BottomHolePressureTarget, TotalRateTarget
    export SinglePhaseRateTarget, DisabledTarget
    export SurfaceLiquidRateTarget, SurfaceOilRateTarget
    export SurfaceWaterRateTarget, SurfaceGasRateTarget
    export HistoricalReservoirVoidageTarget, ReservoirVoidageTarget
    export SinglePhaseRateTarget, BottomHolePressureTarget
    export PerforationMask
    export WellDomain, MultiSegmentWell
    export TotalMassFlux, PotentialDropBalanceWell
    export SegmentWellBoreFrictionHB
    export InjectorControl, ProducerControl
    export Perforations
    export MixedWellSegmentFlow
    export segment_pressure_drop
    export setup_well, setup_vertical_well, setup_well_from_trajectory
    export well_mismatch
    export simulate_data_file, setup_case_from_data_file
    export get_test_setup, get_well_from_mrst_data
    export setup_case_from_mrst
    export simulate_mrst_case
    export MultiPhaseCompositionalSystemLV
    export StandardVolumeSource, VolumeSource, MassSource
    export OverallMoleFractions
    export ImmiscibleSaturation
    export PhaseMassDensities, ConstantCompressibilityDensities
    export BrooksCoreyRelativePermeabilities, TabulatedSimpleRelativePermeabilities
    export EquilibriumRegion
    export setup_reservoir_dict_optimization, optimize_reservoir, parameters_gradient_reservoir

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
    using GeoEnergyIO
    # Artifacts
    using Artifacts
    using LazyArtifacts

    import DataStructures: OrderedDict
    using DocStringExtensions

    timeit_debug_enabled() = Jutul.timeit_debug_enabled()
    function __init__()
        if !haskey(ENV, "JUTULDARCY_PRESERVE_ENV")
            ENV["OPENBLAS_NUM_THREADS"] = 1
            BLAS.set_num_threads(1)
        end
    end

    include("types.jl")
    include("deck_types.jl")
    include("porousmedia_grids.jl")
    include("utils.jl")
    include("state0.jl")
    include("interpolation.jl")
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

    include("flux.jl")
    include("flux_nfvm.jl")

    include("porousmedia.jl")
    # MRST inputs and test cases that use MRST input
    # and .DATA file simulation
    include("input_simulation/input_simulation.jl")
    # Initialization by equilibriation
    include("init/init.jl")
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
    include("coarsening/coarsening.jl")

    include("ext/ext.jl")
    # Nonlinear domain decomposition solvers
    include("NLDD/NLDD.jl")
    # CO2-brine properties
    include("CO2Properties/CO2Properties.jl")
    # Timestepping
    include("timesteps.jl")

    # Geothermal
    include("Geothermal/Geothermal.jl")
    import JutulDarcy.Geothermal: setup_btes_well, setup_vertical_btes_well
    import JutulDarcy.Geothermal: ClosedLoopSupplyToReturnMassCT, ClosedLoopSupplyToReturnEnergyCT, BTESWellGroutEnergyCT
    export setup_btes_well, setup_vertical_btes_well

    # Tracers
    include("Tracers/Tracers.jl")
    import JutulDarcy.Tracers: SinglePhaseTracer, MultiPhaseTracer, add_tracers_to_model!, number_of_tracers
    export SinglePhaseTracer, MultiPhaseTracer, add_tracers_to_model!, number_of_tracers


    @compile_workload begin
        try
            precompile_darcy_multimodels()
            # We run a tiny MRST case to precompile the .MAT file loading
            spe1_path = joinpath(pathof(JutulDarcy), "..", "..", "test", "mrst", "spe1.mat")
            if isfile(spe1_path)
                simulate_mrst_case(spe1_path, info_level = -1, verbose = false)
            end
            # Precompile a DATA file workflow
            spe1_path_data = GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA")
            if isfile(spe1_path_data)
                simulate_data_file(spe1_path_data, info_level = -1)
            end
        catch e
            @warn "Precompilation failure: $e"
        end
    end
end # module
