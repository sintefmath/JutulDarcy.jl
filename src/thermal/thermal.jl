struct BulkVolume <: ScalarVariable end
function Jutul.default_values(model, ::BulkVolume)
    return 1.0
end

function Jutul.default_parameter_values(data_domain, model, param::BulkVolume, symb)
    if haskey(data_domain, :volumes)
        bv = copy(data_domain[:volumes])
    elseif model_or_domain_is_well(model)
        w = physical_representation(data_domain)
        bv = domain_bulk_volume(data_domain, w)
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

"""
    ThermalDispersivity()

Cell-wise thermal dispersivity with two components per cell:
row 1 is longitudinal and row 2 is transversal.
"""
struct ThermalDispersivity <: VectorVariables end
Jutul.values_per_entity(model, ::ThermalDispersivity) = 2
Jutul.minimum_value(::ThermalDispersivity) = 0.0
Jutul.associated_entity(::ThermalDispersivity) = Cells()

function Jutul.default_parameter_values(data_domain, model, param::ThermalDispersivity, symb)
    nc = number_of_cells(data_domain)
    haskey(data_domain, :thermal_dispersivity) || error(":thermal_dispersivity must be present in DataDomain to initialize parameter $symb")
    raw = copy(data_domain[:thermal_dispersivity])

    if raw isa Number
        D = fill(raw, 2, nc)
    elseif raw isa AbstractVector
        n = length(raw)
        if n == nc
            D = zeros(eltype(raw), 2, nc)
            D[1, :] .= raw
            D[2, :] .= raw
        elseif n == 2
            D = zeros(eltype(raw), 2, nc)
            D[1, :] .= raw[1]
            D[2, :] .= raw[2]
        else
            error("Expected :thermal_dispersivity vector length 2 or num_cells ($(nc)), got length $(n)")
        end
    elseif raw isa AbstractMatrix
        size(raw, 1) == 2 || error("Expected :thermal_dispersivity matrix with 2 rows, got size $(size(raw))")
        size(raw, 2) == nc || error("Expected :thermal_dispersivity matrix with size (2, $(nc)), got size $(size(raw))")
        D = raw
    else
        error("Unsupported type $(typeof(raw)) for :thermal_dispersivity")
    end
    return ensure_non_negative_trans(D, "thermal_dispersivity")
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
    fsys = model.system
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
    FluidThermalConductivity()

Cell-wise fluid thermal conductivity.
"""
struct FluidThermalConductivity <: PhaseVariables end
Jutul.default_value(model, ::FluidThermalConductivity) = 0.6
Jutul.minimum_value(::FluidThermalConductivity) = 0.0

function Jutul.default_parameter_values(data_domain, model, param::FluidThermalConductivity, symb)
    nph = number_of_phases(model.system)
    nc = number_of_cells(data_domain)
    if haskey(data_domain, :fluid_thermal_conductivity, Cells())
        C = copy(data_domain[:fluid_thermal_conductivity, Cells()])
    elseif haskey(data_domain, :fluid_thermal_conductivity)
        C = copy(data_domain[:fluid_thermal_conductivity])
    else
        C = default_value(model, param)
    end

    if C isa Number
        λ = fill(C, nph, nc)
    elseif C isa AbstractVector
        if length(C) == nph
            λ = repeat(C, 1, nc)
        else
            length(C) == nc || error("Expected a vector with length num_cells ($(nc)) or num_phases ($(nph)) for :fluid_thermal_conductivity, got length $(length(C)).")
            λ = repeat(C', nph, 1)
        end
    else
        size(C, 1) == nph || error("Expected size $(nph) x num_cells for :fluid_thermal_conductivity, got size $(size(C)).")
        size(C, 2) == nc || error("Expected size $(nph) x $(nc) for :fluid_thermal_conductivity, got size $(size(C)).")
        λ = C
    end
    return ensure_non_negative_trans(λ, "fluid_thermal_conductivity")
end

"""
    FluidThermalTransmissibilites()

Variable defining the fluid component conductivity.
"""
struct FluidThermalTransmissibilites <: VectorVariables end
Jutul.variable_scale(::FluidThermalTransmissibilites) = 1.0
Jutul.minimum_value(::FluidThermalTransmissibilites) = 0.0
Jutul.values_per_entity(model, ::FluidThermalTransmissibilites) = number_of_phases(model.system)

function Jutul.default_parameter_values(data_domain, model, param::FluidThermalTransmissibilites, symb)
    if haskey(data_domain, :fluid_thermal_transmissibilites, Faces())
        # This takes precedence
        T = copy(data_domain[:fluid_thermal_transmissibilites])
    elseif haskey(data_domain, :fluid_thermal_conductivities, Faces())
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
            size(C, 1) == nph || error("Expected size $(nph) x num_cells for :fluid_thermal_conductivity, got size $(size(C))")
            nf = number_of_faces(data_domain)
            T = zeros(nph, nf)
            for ph in 1:nph
                T[ph, :] = compute_face_trans(data_domain, phi.*C[ph, :])
            end
        end
    elseif haskey(data_domain, :fluid_thermal_conductivity)
        nph = number_of_phases(model.system)
        C = data_domain[:fluid_thermal_conductivity]
        phi = data_domain[:porosity]
        if C isa Vector
            T = compute_face_trans(data_domain, phi.*C)
            T = repeat(T', nph, 1)
        else
            size(C, 1) == nph || error("Expected size $(nph) x num_cells for :fluid_thermal_conductivity, got size $(size(C))")
            nf = number_of_faces(data_domain)
            T = zeros(nph, nf)
            for ph in 1:nph
                T[ph, :] = compute_face_trans(data_domain, phi.*C[ph, :])
            end
        end
    else
        nph = number_of_phases(model.system)
        C = default_value(model, FluidThermalConductivity())
        phi = data_domain[:porosity]
        T = compute_face_trans(data_domain, phi.*fill(C, number_of_cells(data_domain)))
        T = repeat(T', nph, 1)
    end
    return ensure_non_negative_trans(T, "fluid_thermal_transmissibilites")
end

Jutul.associated_entity(::FluidThermalTransmissibilites) = Faces()

struct RockThermalTransmissibilites <: ScalarVariable end
Jutul.variable_scale(::RockThermalTransmissibilites) = 1.0
Jutul.minimum_value(::RockThermalTransmissibilites) = 0.0
Jutul.associated_entity(::RockThermalTransmissibilites) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::RockThermalTransmissibilites, symb)
    if haskey(data_domain, :rock_thermal_transmissibilites, Faces())
        # This takes precedence
        T = copy(data_domain[:rock_thermal_transmissibilites])
    elseif haskey(data_domain, :rock_thermal_conductivities, Faces())
        # This takes precedence
        T = copy(data_domain[:rock_thermal_conductivities])
    elseif haskey(data_domain, :rock_thermal_conductivity, Cells())
        T = reservoir_conductivity(data_domain)
    else
        error(":rock_thermal_transmissibilites/:rock_thermal_conductivities or :rock_thermal_conductivity symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return ensure_non_negative_trans(T, "rock_thermal_transmissibilites")
end

const FluidThermalConductivities = FluidThermalTransmissibilites
const RockThermalConductivities = RockThermalTransmissibilites

function ensure_non_negative_trans(T, name)
    bad = 0
    neg = 0
    for (i, v) in enumerate(T)
        if !isfinite(v)
            T[i] = 0.0
            bad += 1
        elseif v < 0.0
            T[i] = 0.0
            neg += 1
        end
    end
    if neg > 0
        jutul_message(name, "Found $neg negative values, set to zero.")
    end
    if bad > 0
        jutul_message(name, "Found $bad non-finite values, set to zero.")
    end
    return T
end

function reservoir_conductivity(reservoir::DataDomain)
    phi = reservoir[:porosity]
    C = reservoir[:rock_thermal_conductivity]
    T = compute_face_trans(reservoir, (1.0 .- phi).*C)
    if haskey(reservoir, :nnc)
        nnc = reservoir[:nnc]
        nnc::NonNeighboringConnections
        num_nnc = length(nnc.trans_thermal)
        # NNC come at the end.
        offset = number_of_faces(reservoir) - num_nnc
        for (i, T_nnc) in enumerate(nnc.trans_thermal)
            T[i + offset] = T_nnc
        end
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

function Jutul.default_parameter_values(data_domain, model, param::WellIndicesThermal, symb)

    WIt = copy(data_domain[:thermal_well_index, Perforations()])
    dims = data_domain[:cell_dims, Perforations()]
    thermal_conductivity = data_domain[:thermal_conductivity, Perforations()]
    direction = data_domain[:perforation_direction, Perforations()]
    radius = data_domain[:perforation_radius, Perforations()]
    drainage_radius = data_domain[:drainage_radius, Perforations()]
    gdim = size(data_domain[:cell_centroids, Cells()], 1
)
    # These are defined per cell, map to perforations
    well = physical_representation(data_domain)
    ic = well.perforations.self

    λ_casing = data_domain[:casing_thermal_conductivity, Cells()][ic]
    λ_grout = data_domain[:grouting_thermal_conductivity, Cells()][ic]
    casing_thickness = data_domain[:casing_thickness, Cells()][ic]
    grouting_thickness = data_domain[:grouting_thickness, Cells()][ic]

    T = Base.promote_type(
        eltype(WIt), eltype(thermal_conductivity), eltype(radius),
        eltype(λ_casing), eltype(λ_grout), eltype(casing_thickness),
        eltype(grouting_thickness), eltype(drainage_radius))
    if T != eltype(WIt)
        WIt = convert(Vector{T}, WIt)
    end

    for (i, val) in enumerate(WIt)
        defaulted = !isfinite(val)
        if defaulted
            Δ = dims[i]
            if thermal_conductivity isa AbstractVector
                Λ_i = thermal_conductivity[i]
            else
                Λ_i = thermal_conductivity[:, i]
            end
            Λ_i = Jutul.expand_perm(Λ_i, gdim)
            WIt[i] = compute_well_thermal_index(Δ, Λ_i, radius[i], direction[i];
                casing_thickness = casing_thickness[i],
                grouting_thickness = grouting_thickness[i],
                casing_thermal_conductivity = λ_casing[i],
                grouting_thermal_conductivity = λ_grout[i],
                drainage_radius = drainage_radius[i],
            )
        end
    end
    return WIt
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
    MaterialDensities()

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
        error(":material_density must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
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
        error(":material_heat_capacity symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
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
            RockThermalTransmissibilites = RockThermalTransmissibilites(),
            FluidThermalTransmissibilites = FluidThermalTransmissibilites(),
            FluidThermalConductivity = FluidThermalConductivity()
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
    set_primary_variables!(facility,
        SurfaceTemperature = SurfaceTemperature(),
        SurfaceEnthalpy = SurfaceEnthalpy()
    )
    facility.equations[:temperature_equation] = SurfaceTemperatureEquation()
    facility.equations[:enthalpy_equation] = SurfaceEnthalpyEquation()
    out = facility.output_variables
    push!(out, :SurfaceTemperature)
    push!(out, :SurfaceEnthalpy)
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

function model_is_thermal(model::SimulationModel, return_name=false)
    pvars = Jutul.get_primary_variables(model)
    candidates = [:Temperature, :Enthalpy]
    is_thermal, name = false, nothing
    for var in candidates
        if haskey(pvars, var)
            is_thermal = true
            name = var
            break
        end
    end
    if return_name
        return is_thermal, name
    else
        return is_thermal
    end
end

include("variables.jl")
include("equations.jl")