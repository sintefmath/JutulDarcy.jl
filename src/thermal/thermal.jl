export ThermalSystem
struct ThermalSystem{T} <: JutulSystem
    nph::Int64
    function ThermalSystem(; nphases = 1, formulation = :Temperature)
        @assert formulation == :Temperature
        new{formulation}(nphases)
    end
end

thermal_system(sys::ThermalSystem) = sys
thermal_system(sys::CompositeSystem) = sys.systems.thermal

const ThermalModel = SimulationModel{<:Any, <:ThermalSystem, <:Any, <:Any}

struct BulkVolume <: ScalarVariable end
function Jutul.default_values(model, ::BulkVolume)
    return 5*fluid_volume(model.domain)
end

struct RockHeatCapacity <: ScalarVariable end
Jutul.default_value(model, ::RockHeatCapacity) = 5000.0
struct RockDensity <: ScalarVariable end
Jutul.default_value(model, ::RockDensity) = 1500.0

struct RockInternalEnergy <: ScalarVariable end
struct TotalThermalEnergy <: ScalarVariable end

struct FluidHeatCapacity <: PhaseVariables end
Jutul.default_value(model, ::FluidHeatCapacity) = 10000.0
struct FluidInternalEnergy <: PhaseVariables end
struct FluidEnthalpy <: PhaseVariables end

struct FluidThermalConductivities <: ScalarVariable end
Jutul.variable_scale(::FluidThermalConductivities) = 1e-10
Jutul.minimum_value(::FluidThermalConductivities) = 0.0
Jutul.default_value(model, ::FluidThermalConductivities) = 1e-3
Jutul.associated_entity(::FluidThermalConductivities) = Faces()

struct RockThermalConductivities <: ScalarVariable end
Jutul.variable_scale(::RockThermalConductivities) = 1e-10
Jutul.minimum_value(::RockThermalConductivities) = 0.0
Jutul.default_value(model, ::RockThermalConductivities) = 1e-3
Jutul.associated_entity(::RockThermalConductivities) = Faces()

number_of_phases(t::ThermalSystem) = t.nph

function select_primary_variables!(S, system::ThermalSystem, model)
    S[:Temperature] = Temperature()
end

function select_parameters!(S, system::ThermalSystem, model)
    nph = number_of_phases(system)
    # Rock itself
    S[:RockHeatCapacity] = RockHeatCapacity()
    S[:RockDensity] = RockDensity()
    S[:BulkVolume] = BulkVolume()
    # Fluid heat related parameters
    S[:FluidHeatCapacity] = FluidHeatCapacity()
    S[:FluidVolume] = FluidVolume()
    # Fluid flow related parameters
    S[:PhaseMassDensities] = ConstantCompressibilityDensities(nph)
    S[:Pressure] = Pressure()
    S[:PhaseMassMobilities] = PhaseMassMobilities()
end

function select_parameters!(prm, disc::D, model::ThermalModel) where D<:Union{TwoPointPotentialFlowHardCoded, Jutul.PotentialFlow}
    prm[:FluidThermalConductivities] = FluidThermalConductivities()
    prm[:RockThermalConductivities] = RockThermalConductivities()
    prm[:Transmissibilities] = Transmissibilities()
    prm[:TwoPointGravityDifference] = TwoPointGravityDifference()
end

function select_secondary_variables!(S, system::ThermalSystem, model)
    nph = number_of_phases(system)
    S[:FluidInternalEnergy] = FluidInternalEnergy()
    S[:FluidEnthalpy] = FluidEnthalpy()
    S[:RockInternalEnergy] = RockInternalEnergy()
    S[:TotalThermalEnergy] = TotalThermalEnergy()
end

include("variables.jl")
include("equations.jl")