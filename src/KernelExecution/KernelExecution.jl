module KernelExecution
    import Adapt
    using Jutul
    import Jutul.KernelExecution: KernelAbstractions
    using .KernelAbstractions: @Const, @index, @kernel
    using LinearAlgebra
    using SparseArrays: nonzeros
    using StaticArrays

    import ..JutulDarcy: DeckPhaseMassDensities, DeckPhaseViscosities,
        DeckShrinkageFactors, DeckThermalViscosityTable, FacilitySystem,
        LinearlyCompressiblePoreVolume, MinimalTPFATopology, MuBTable,
        PerforationMask, PhaseRelativePermeability, PVDG, PVDO, PVTG,
        PVTGTable, PVTO, PVTOTable, PVTW, ReservoirFromWellFlowCT,
        ReservoirFromWellThermalCT, ReservoirRelativePermeabilities, Rs, Rv,
        ScaledCapillaryPressure, SimpleCapillaryPressure, SimpleWell,
        StandardBlackOilSystem, TemperatureDependentVariable, WellGroup,
        cpr_weights_no_partials!, kernel_abstractions_backend,
        update_analytical_cpr_weights!, update_p_rhs!,
        update_pressure_system!, update_quasi_impes_weights!,
        update_true_impes_weights!

    kernel_abstractions_backend(::Val{:ka}) = KernelAbstractions.CPU()
    kernel_abstractions_backend(::Val{:ka_cpu}) = KernelAbstractions.CPU()

    include("adapt.jl")
    include("cpr.jl")
end
