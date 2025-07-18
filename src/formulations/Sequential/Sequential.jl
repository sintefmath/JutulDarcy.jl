module Sequential
    import Jutul
    import JutulDarcy

    import JutulDarcy:
        MultiPhaseSystem,
        SourceTerm,
        FlowBoundaryCondition,
        SimpleWellModel,
        PhaseVariables,
        Pressure,
        Saturations,
        ImmiscibleModel,
        StandardBlackOilModel,
        StandardBlackOilModelWithWater,
        ReservoirFromWellFlowCT,
        BlackOilVariableSwitchingSystem,
        Rs,
        Rv,
        has_vapoil,
        has_disgas,
        number_of_phases,
        number_of_components,
        reference_densities,
        reservoir_linsolve,
        default_psolve,
        phase_source,
        reservoir_model,
        compute_bc_mass_fluxes,
        apply_flow_bc!,
        darcy_phase_volume_fluxes,
        cross_term_perforation_get_conn,
        multisegment_well_perforation_flux!,
        effective_transmissibility,
        effective_gravity_difference,
        face_average_density,
        capillary_pressure,
        capillary_gradient,
        model_or_domain_is_well

    import Jutul:
        JutulCase,
        MultiModel,
        SimulationModel,
        ConservationLaw,
        JutulStorage,
        JutulEquation,
        Simulator,
        JutulSimulator,
        LocalStateAD,
        JutulFormulation,
        AdditiveCrossTerm,
        FluxType,
        ScalarVariable,
        TwoPointPotentialFlowHardCoded,
        EquationMajorLayout,
        ParallelCSRContext,
        Cells,
        Faces,
        TPFA,
        SPU,
        CompactAutoDiffCache,
        ProgressRecorder,
        GenericKrylov,
        @tic,
        @jutul_secondary,
        number_of_cells,
        number_of_faces,
        number_of_half_faces,
        degrees_of_freedom_per_entity,
        number_of_equations_per_entity,
        number_of_equations,
        number_of_entities,
        count_active_entities,
        global_map,
        update_secondary_variable!,
        get_dependencies,
        conserved_symbol,
        get_entries,
        setup_parameters,
        setup_state,
        setup_state_and_parameters,
        local_discretization,
        value,
        as_value,
        cell_pair,
        update_values!,
        jutul_message

    import StaticArrays: SVector

    include("types.jl")
    include("utils.jl")
    include("pressure/pressure.jl")
    include("transport/transport.jl")
end
