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

function Jutul.subvariable(p::TemperatureDependentVariable, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(p.regions, c)
    return TemperatureDependentVariable(p.tab, regions = regions)
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

function Jutul.subvariable(p::PressureTemperatureDependentVariable, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(p.regions, c)
    return PressureTemperatureDependentVariable(p.tab, regions = regions)
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

@jutul_secondary function update_temperature_dependent_enthalpy!(H_phases, var::PressureTemperatureDependentEnthalpy{T, R, N}, model::CompositionalModel, Pressure, Temperature, LiquidMassFractions, VaporMassFractions, PhaseMassDensities, ix) where {T, R, N}
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
        phi = data_domain[:porosity]
        if C isa Vector
            T = compute_face_trans(data_domain, phi.*C)
            T = repeat(T', nph, 1)
        else
            @assert size(C, 1) == nph
            nf = number_of_faces(data_domain)
            T = zeros(nph, nf)
            for ph in 1:nph
                T[ph, :] = compute_face_trans(data_domain, phi.*C[ph, :])
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
        phi = data_domain[:porosity]
        C = data_domain[:rock_thermal_conductivity]
        T = compute_face_trans(data_domain, (1.0 .- phi).*C)
    else
        error(":rock_thermal_conductivities or :rock_thermal_conductivities symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end


"""
    WellIndicesThermal()

Parameter for the thermal connection strength between a well and the reservoir
for a given perforation. Typical values come from a combination of Peaceman's
formula with thermal conducivity in place of permeability, upscaling and/or
history matching.
"""
struct WellIndicesThermal <: ScalarVariable end

Jutul.minimum_value(::WellIndicesThermal) = 0.0
Jutul.variable_scale(::WellIndicesThermal) = 1.0
Jutul.associated_entity(::WellIndicesThermal) = Perforations()
function Jutul.default_values(model, ::WellIndicesThermal)
    w = physical_representation(model.domain)
    return vec(copy(w.perforations.WIth))
end

"""
    MaterialThermalConductivities()

Parameter for the thermal conductivities of the materials in the well.
"""
struct MaterialThermalConductivities <: ScalarVariable end

Jutul.variable_scale(::MaterialThermalConductivities) = 1e-10
Jutul.minimum_value(::MaterialThermalConductivities) = 0.0
Jutul.associated_entity(::MaterialThermalConductivities) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::MaterialThermalConductivities, symb)
    if haskey(data_domain, :material_thermal_conductivity, Faces())
        T = copy(data_domain[:material_thermal_conductivity])
    else
        error(":material_thermal_conductivity or :material_thermal_conductivity symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

"""
    MaterialDensity()

Parameter well material density.
"""
struct MaterialDensities <: ScalarVariable end

Jutul.variable_scale(::MaterialDensities) = 1.0
Jutul.minimum_value(::MaterialDensities) = 0.0
Jutul.associated_entity(::MaterialDensities) = Cells()

function Jutul.default_parameter_values(data_domain, model, param::MaterialDensities, symb)
    if haskey(data_domain, :material_density, Cells())
        T = copy(data_domain[:material_density])
    else
        error(":material_density or :material_density symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

"""
    MaterialHeatCapacities()

Parameter heat capacitiy of the well material.
"""
struct MaterialHeatCapacities <: ScalarVariable end

Jutul.variable_scale(::MaterialHeatCapacities) = 1.0
Jutul.minimum_value(::MaterialHeatCapacities) = 0.0
Jutul.associated_entity(::MaterialHeatCapacities) = Cells()

function Jutul.default_parameter_values(data_domain, model, param::MaterialHeatCapacities, symb)
    if haskey(data_domain, :material_heat_capacity, Cells())
        T = copy(data_domain[:material_heat_capacity])
    else
        error(":material_heat_capacity or :material_heat_capacity symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end


struct MaterialInternalEnergy <: ScalarVariable end

"""
    add_thermal_to_model!(model::MultiModel)

Add energy conservation equation and thermal primary variable together with
standard set of parameters to existing flow model. Note that more complex models
require additional customization after this function call to get correct
results.
"""
function add_thermal_to_model!(model::MultiModel)
    for (k, m) in pairs(model.models)
        if m.system isa MultiPhaseSystem
            add_thermal_to_model!(m)
        elseif m.system isa FacilitySystem
            add_thermal_to_facility!(m)
        end
    end
    return m
end

function add_thermal_to_model!(model)
    set_primary_variables!(model, Temperature = Temperature())
    set_parameters!(model,
        RockHeatCapacity = RockHeatCapacity(),
        RockDensity = RockDensity(),
        BulkVolume = BulkVolume(),
        ComponentHeatCapacity = ComponentHeatCapacity(),
    )
    set_secondary_variables!(model,
        FluidInternalEnergy = FluidInternalEnergy(),
        FluidEnthalpy = FluidEnthalpy(),
        TotalThermalEnergy = TotalThermalEnergy(),
    )
    is_reservoir = !model_or_domain_is_well(model)
    if is_reservoir
        set_parameters!(model,
            RockThermalConductivities = RockThermalConductivities(),
            FluidThermalConductivities = FluidThermalConductivities()
        )
        set_secondary_variables!(model,
            RockInternalEnergy = RockInternalEnergy()
        )
    else
        if model_or_domain_is_well(model)
            w = physical_representation(model.domain)

            set_parameters!(model,
                WellIndicesThermal = WellIndicesThermal(),
            )
            if w isa MultiSegmentWell
                set_parameters!(model,
                    MaterialThermalConductivities = MaterialThermalConductivities(),
                    MaterialHeatCapacities = MaterialHeatCapacities(),
                    MaterialDensities = MaterialDensities()
                )
                set_secondary_variables!(model,
                    MaterialInternalEnergy = MaterialInternalEnergy()
                )

            else
                w::SimpleWell
                set_secondary_variables!(model,
                    RockInternalEnergy = RockInternalEnergy()
                )
            end
        end
    end
    disc = model.domain.discretizations.heat_flow
    model.equations[:energy_conservation] = ConservationLaw(disc, :TotalThermalEnergy, 1)

    out = model.output_variables
    push!(out, :TotalThermalEnergy)
    push!(out, :FluidEnthalpy)
    push!(out, :Temperature)

    unique!(out)
    return model
end

function add_thermal_to_facility!(facility)
    set_primary_variables!(facility, SurfaceTemperature = SurfaceTemperature())
    facility.equations[:temperature_equation] = SurfaceTemperatureEquation()
    out = facility.output_variables
    push!(out, :SurfaceTemperature)
    unique!(out)
    return facility
end

"""
    model_is_thermal(model)

Utility function to check if a model has thermal equations.
"""
function model_is_thermal(model::MultiModel)
    m = reservoir_model(model)
    return model_is_thermal(m)
end

function model_is_thermal(model::SimulationModel)
    pvars = Jutul.get_primary_variables(model)
    return haskey(pvars, :Temperature)
end

include("variables.jl")
include("equations.jl")