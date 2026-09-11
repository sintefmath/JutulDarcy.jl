module KernelExecution
    import Adapt
    using Jutul
    import Jutul.KernelExecution: KernelAbstractions

    import ..JutulDarcy: DeckPhaseMassDensities, DeckPhaseViscosities,
        DeckShrinkageFactors, DeckThermalViscosityTable, FacilitySystem,
        LinearlyCompressiblePoreVolume, MinimalTPFATopology, MuBTable,
        PerforationMask, PhaseRelativePermeability, PVDG, PVDO, PVTG,
        PVTGTable, PVTO, PVTOTable, PVTW, ReservoirFromWellFlowCT,
        ReservoirFromWellThermalCT, ReservoirRelativePermeabilities, Rs, Rv,
        ScaledCapillaryPressure, SimpleCapillaryPressure, SimpleWell,
        StandardBlackOilSystem, TemperatureDependentVariable, WellGroup,
        kernel_abstractions_backend

    kernel_abstractions_backend(::Val{:ka}) = KernelAbstractions.CPU()
    kernel_abstractions_backend(::Val{:ka_cpu}) = KernelAbstractions.CPU()

    include("adapt.jl")
end
