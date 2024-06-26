

"""
    ThermalSystem(flow_system, formulation = :Temperature)

Geothermal system that defines heat transfer through fluid advection and through
the rock itself. Can be combined with a multiphase system using
[`Jutul.CompositeSystem`](@ref).
"""
struct ThermalSystem{T, F} <: JutulSystem
    flow_system::F
    function ThermalSystem(sys::T;
            formulation = :Temperature
        ) where T<:Union{MultiPhaseSystem, Missing}
        @assert formulation == :Temperature
        new{formulation, T}(sys)
    end
end

thermal_system(sys::ThermalSystem) = sys
thermal_system(sys::CompositeSystem) = sys.systems.thermal
flow_system(sys::ThermalSystem) = sys.flow_system

const ThermalModel = SimulationModel{<:JutulDomain, <:ThermalSystem, <:Any, <:Any}

const ThermalCompositionalModel = SimulationModel{<:JutulDomain, <:ThermalSystem{<:Any, <:CompositionalSystem}, <:Any, <:Any}
const ThermalBlackOilModel = SimulationModel{<:JutulDomain, <:ThermalSystem{<:Any, <:StandardBlackOilModel}, <:Any, <:Any}
const ThermalImmiscibleModel = SimulationModel{<:JutulDomain, <:ThermalSystem{<:Any, <:ImmiscibleSystem}, <:Any, <:Any}
const ThermalSinglePhaseModel = SimulationModel{<:JutulDomain, <:ThermalSystem{<:Any, <:SinglePhaseSystem}, <:Any, <:Any}

struct BulkVolume <: ScalarVariable end
function Jutul.default_values(model, ::BulkVolume)
    return 1.0
end

function Jutul.default_parameter_values(data_domain, model, param::BulkVolume, symb)
    if haskey(data_domain, :volumes)
        bv = copy(data_domain[:volumes])
    elseif model_or_domain_is_well(model)
        r = physical_representation(model.domain)
        if r isa MultiSegmentWell
            bv = copy(model.domain.representation.volumes)
        else
            r::SimpleWell
            bv = [r.volume]
        end
    end
    return bv
end

struct RockHeatCapacity <: ScalarVariable end
Jutul.default_value(model, ::RockHeatCapacity) = 1000.0

function Jutul.default_parameter_values(data_domain, model, param::RockHeatCapacity, symb)
    if haskey(data_domain, :rock_heat_capacity, Cells())
        # This takes precedence
        T = copy(data_domain[:rock_heat_capacity])
    else
        T = fill(default_value(model, param), number_of_cells(data_domain))
    end
    return T
end

struct RockDensity <: ScalarVariable end
Jutul.default_value(model, ::RockDensity) = 2000.0

function Jutul.default_parameter_values(data_domain, model, param::RockDensity, symb)
    if haskey(data_domain, :rock_density, Cells())
        # This takes precedence
        T = copy(data_domain[:rock_density])
    else
        T = fill(default_value(model, param), number_of_cells(data_domain))
    end
    return T
end

struct RockInternalEnergy <: ScalarVariable end
struct TotalThermalEnergy <: ScalarVariable end

struct ComponentHeatCapacity <: ComponentVariables end
Jutul.default_value(model, ::ComponentHeatCapacity) = 4184.0

function Jutul.default_parameter_values(data_domain, model, param::ComponentHeatCapacity, symb)
    ncomp = number_of_components(model.system)
    if haskey(data_domain, :component_heat_capacity, Cells())
        # This takes precedence
        T = copy(data_domain[:component_heat_capacity])
        if T isa Vector
            T = repeat(T', ncomp, 1)
        else
            @assert size(T, 1) == ncomp
        end
    else
        T = fill(default_value(model, param), ncomp, number_of_cells(data_domain))
    end
    return T
end

struct FluidInternalEnergy <: PhaseVariables end
struct FluidEnthalpy <: PhaseVariables end

struct TemperatureDependentVariable{T, R, N} <: VectorVariables
    tab::T
    regions::R
    function TemperatureDependentVariable(tab; regions = nothing)
        tab = region_wrap(tab, regions)
        ex = first(tab)
        N = length(ex(273.15 + 30.0))
        new{typeof(tab), typeof(regions), N}(tab, regions)
    end
end

function Jutul.values_per_entity(model, ::TemperatureDependentVariable{T, R, N}) where {T, R, N}
    return N
end

@jutul_secondary function update_temperature_dependent!(result, var::TemperatureDependentVariable{T, R, N}, model, Temperature, ix) where {T, R, N}
    for c in ix
        reg = region(var.regions, c)
        interpolator = table_by_region(var.tab, reg)
        F_of_T = interpolator(Temperature[c])
        for i in 1:N
            result[i, c] = F_of_T[i]
        end
    end
    return result
end

struct PressureTemperatureDependentVariable{T, R, N} <: VectorVariables
    tab::T
    regions::R
    function PressureTemperatureDependentVariable(tab; regions = nothing)
        tab = region_wrap(tab, regions)
        ex = first(tab)
        N = length(ex(1e8, 273.15 + 30.0))
        new{typeof(tab), typeof(regions), N}(tab, regions)
    end
end

function Jutul.values_per_entity(model, ::PressureTemperatureDependentVariable{T, R, N}) where {T, R, N}
    return N
end

@jutul_secondary function update_temperature_dependent!(result, var::PressureTemperatureDependentVariable{T, R, N}, model, Pressure, Temperature, ix) where {T, R, N}
    for c in ix
        reg = region(var.regions, c)
        interpolator = table_by_region(var.tab, reg)
        F_of_T = interpolator(Pressure[c], Temperature[c])
        for i in 1:N
            result[i, c] = F_of_T[i]
        end
    end
    return result
end

struct PressureTemperatureDependentEnthalpy{T, R, N} <: VectorVariables
    tab::T
    regions::R
    function PressureTemperatureDependentEnthalpy(tab; regions = nothing)
        tab = region_wrap(tab, regions)
        ex = first(tab)
        N = length(ex(1e8, 273.15 + 30.0))
        new{typeof(tab), typeof(regions), N}(tab, regions)
    end
end

function Jutul.values_per_entity(model, ::PressureTemperatureDependentEnthalpy{T, R, N}) where {T, R, N}
    return N
end

@jutul_secondary function update_temperature_dependent_enthalpy!(H_phases, var::PressureTemperatureDependentEnthalpy{T, R, N}, model::ThermalCompositionalModel, Pressure, Temperature, LiquidMassFractions, VaporMassFractions, PhaseMassDensities, ix) where {T, R, N}
    fsys = flow_system(model.system)
    @assert !has_other_phase(fsys)
    @assert N == number_of_components(fsys)
    l, v = phase_indices(fsys)

    X, Y = LiquidMassFractions, VaporMassFractions
    rho = PhaseMassDensities
    for c in ix
        reg = region(var.regions, c)
        interpolator = table_by_region(var.tab, reg)
        component_H = interpolator(Pressure[c], Temperature[c])
        H_l = 0.0
        H_v = 0.0
        for i in 1:N
            H_i = component_H[i]
            H_l += X[i, c]*H_i
            H_v += Y[i, c]*H_i
        end
        # p = Pressure[c]
        H_phases[l, c] = H_l# + p/rho[l, c]
        H_phases[v, c] = H_v# + p/rho[v, c]
    end
    return H_phases
end

"""
    FluidThermalConductivities()

Variable defining the fluid component conductivity.
"""
struct FluidThermalConductivities <: VectorVariables end
Jutul.variable_scale(::FluidThermalConductivities) = 1e-10
Jutul.minimum_value(::FluidThermalConductivities) = 0.0
Jutul.values_per_entity(model, ::FluidThermalConductivities) = number_of_phases(model.system)

function Jutul.default_parameter_values(data_domain, model, param::FluidThermalConductivities, symb)
    if haskey(data_domain, :fluid_thermal_conductivities, Faces())
        # This takes precedence
        T = copy(data_domain[:fluid_thermal_conductivities])
    elseif haskey(data_domain, :fluid_thermal_conductivity, Cells())
        nph = number_of_phases(model.system)
        C = data_domain[:fluid_thermal_conductivity]
        if C isa Vector
            T = compute_face_trans(data_domain, C)
            T = repeat(T', nph, 1)
        else
            @assert size(C, 1) == nph
            nf = number_of_faces(data_domain)
            T = zeros(nph, nf)
            for ph in 1:nph
                T[ph, :] = compute_face_trans(data_domain, C[ph, :])
            end
        end
    else
        error(":fluid_thermal_conductivities or :fluid_thermal_conductivities symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

Jutul.associated_entity(::FluidThermalConductivities) = Faces()

struct RockThermalConductivities <: ScalarVariable end
Jutul.variable_scale(::RockThermalConductivities) = 1e-10
Jutul.minimum_value(::RockThermalConductivities) = 0.0
Jutul.associated_entity(::RockThermalConductivities) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::RockThermalConductivities, symb)
    if haskey(data_domain, :rock_thermal_conductivities, Faces())
        # This takes precedence
        T = copy(data_domain[:rock_thermal_conductivities])
    elseif haskey(data_domain, :rock_thermal_conductivity, Cells())
        nph = number_of_phases(model.system)
        C = data_domain[:rock_thermal_conductivity]
        T = compute_face_trans(data_domain, C)
    else
        error(":fluid_thermal_conductivities or :fluid_thermal_conductivities symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

number_of_phases(t::ThermalSystem) = number_of_phases(flow_system(t))
number_of_components(t::ThermalSystem) = number_of_components(flow_system(t))

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
    S[:ComponentHeatCapacity] = ComponentHeatCapacity()
    S[:FluidVolume] = FluidVolume()
    # Fluid flow related parameters
    S[:PhaseMassDensities] = ConstantCompressibilityDensities(nph)
    S[:Pressure] = Pressure()
    S[:Saturations] = Saturations()
    S[:PhaseViscosities] = PhaseViscosities()
    if !model_or_domain_is_well(model)
        S[:PhaseMassMobilities] = PhaseMassMobilities()
        S[:RelativePermeabilities] = RelativePermeabilitiesParameter()
    end
end

function select_parameters!(prm, disc::D, model::ThermalModel) where D<:Union{TwoPointPotentialFlowHardCoded, Jutul.PotentialFlow}
    if !model_or_domain_is_well(model)
        prm[:FluidThermalConductivities] = FluidThermalConductivities()
        prm[:RockThermalConductivities] = RockThermalConductivities()
        prm[:Transmissibilities] = Transmissibilities()
        prm[:TwoPointGravityDifference] = TwoPointGravityDifference()
    end
end

function select_secondary_variables!(S, system::ThermalSystem, model)
    S[:FluidInternalEnergy] = FluidInternalEnergy()
    S[:FluidEnthalpy] = FluidEnthalpy()
    S[:RockInternalEnergy] = RockInternalEnergy()
    S[:TotalThermalEnergy] = TotalThermalEnergy()
end

function select_minimum_output_variables!(out, system::ThermalSystem, model)
    push!(out, :TotalThermalEnergy)
    push!(out, :FluidEnthalpy)
end

include("variables.jl")
include("equations.jl")