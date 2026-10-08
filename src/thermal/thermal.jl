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

struct RockInternalEnergy <: ScalarVariable end
struct TotalThermalEnergy <: ScalarVariable end

"""
    PotentialEnergy()

Gravitational potential energy of the fluid in each cell, i.e. the total fluid
mass times [`UnitPotentialEnergy`](@ref). Only used with the `:total` energy
formulation.
"""
struct PotentialEnergy <: ScalarVariable end

"""
    TotalEnergy()

Total energy in each cell, defined as the sum of `TotalThermalEnergy`
and [`PotentialEnergy`](@ref). This is the conserved quantity when the `:total`
energy formulation is used.
"""
struct TotalEnergy <: ScalarVariable end

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

"""
    FluidInternalEnergy()

Specific internal energy of each fluid phase (J/kg). The default implementation
is `C*T` where `C` is `ComponentHeatCapacity` and `T` the temperature.
See [`PressureTemperatureDependentInternalEnergy`](@ref) for a tabulated
alternative.
"""
struct FluidInternalEnergy <: PhaseVariables end
struct FluidEnthalpy <: PhaseVariables end

"""
    PressureTemperatureDependentInternalEnergy(tab; regions = nothing)

Specific internal energy of each fluid phase (J/kg) given by a table `tab` that
is evaluated as `tab(p, T)`. Used in place of [`FluidInternalEnergy`](@ref),
replacing the default `C*T` implementation.

For immiscible and single-phase systems, the table gives one value per phase.
For compositional systems, the table gives one value per component (pure
component internal energy), and the phase values are obtained by mass fraction
weighting (ideal mixing). If the system has an additional immiscible aqueous
phase, the first entry in the table corresponds to that phase.

When this variable is used, the default enthalpy of injected fluid (see
[`InjectorControl`](@ref)) is evaluated from the same table at the injection
temperature, ensuring that the injected enthalpy uses the same reference state
as the rest of the model.
"""
struct PressureTemperatureDependentInternalEnergy{T, R} <: PhaseVariables
    tab::T
    regions::R
    function PressureTemperatureDependentInternalEnergy(tab; regions = nothing)
        tab = region_wrap(tab, regions)
        new{typeof(tab), typeof(regions)}(tab, regions)
    end
    function PressureTemperatureDependentInternalEnergy(tab::T, regions::R,
            ::Val{:assembled}) where {T, R}
        # Internal constructor for already processed tables (e.g. for Adapt)
        return new{T, R}(tab, regions)
    end
end

function Jutul.subvariable(p::PressureTemperatureDependentInternalEnergy, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(p.regions, c)
    return PressureTemperatureDependentInternalEnergy(p.tab, regions = regions)
end

"""
    tabulated_internal_energy(var::PressureTemperatureDependentInternalEnergy, p, T, cell)

Evaluate the table of `var` at pressure `p` and temperature `T` for the region
of `cell`. Returns per-phase values (immiscible) or per-component values
(compositional).
"""
function tabulated_internal_energy(var::PressureTemperatureDependentInternalEnergy, p, T, cell)
    interpolator = table_by_region(var.tab, region(var.regions, cell))
    return interpolator(p, T)
end

struct TemperatureDependentVariable{T, R, N} <: VectorVariables
    tab::T
    regions::R
    function TemperatureDependentVariable(tab; regions = nothing)
        tab = region_wrap(tab, regions)
        ex = first(tab)
        N = length(ex(273.15 + 30.0))
        new{typeof(tab), typeof(regions), N}(tab, regions)
    end
    function TemperatureDependentVariable(tab::T, regions::R,
            ::Val{N}) where {T, R, N}
        new{T, R, N}(tab, regions)
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
    function PressureTemperatureDependentVariable(tab::T, regions::R,
            ::Val{N}) where {T, R, N}
        return new{T, R, N}(tab, regions)
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
    function PressureTemperatureDependentEnthalpy(tab::T, regions::R,
            ::Val{N}) where {T, R, N}
        return new{T, R, N}(tab, regions)
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
    FluidThermalConductivities()

Variable defining the fluid component conductivity.
"""
struct FluidThermalConductivities <: VectorVariables end
Jutul.variable_scale(::FluidThermalConductivities) = 1.0
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
            size(C, 1) == nph || error("Expected size $(nph) x num_cells for :fluid_thermal_conductivity, got size $(size(C))")
            nf = number_of_faces(data_domain)
            T = zeros(nph, nf)
            for ph in 1:nph
                T[ph, :] = compute_face_trans(data_domain, phi.*C[ph, :])
            end
        end
    else
        error(":fluid_thermal_conductivities or :fluid_thermal_conductivity symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return ensure_non_negative_trans(T, "fluid_thermal_conductivities")
end

Jutul.associated_entity(::FluidThermalConductivities) = Faces()

struct RockThermalConductivities <: ScalarVariable end
Jutul.variable_scale(::RockThermalConductivities) = 1.0
Jutul.minimum_value(::RockThermalConductivities) = 0.0
Jutul.associated_entity(::RockThermalConductivities) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::RockThermalConductivities, symb)
    if haskey(data_domain, :rock_thermal_conductivities, Faces())
        # This takes precedence
        T = copy(data_domain[:rock_thermal_conductivities])
    elseif haskey(data_domain, :rock_thermal_conductivity, Cells())
        T = reservoir_conductivity(data_domain)
    else
        error(":rock_thermal_conductivities or :rock_thermal_conductivity symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return ensure_non_negative_trans(T, "rock_thermal_conductivities")
end

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
    UnitPotentialEnergy()

Parameter for the gravitational potential energy per unit mass in each cell
(J/kg). The default is `-g*z` where `z` is the depth of the cell centroid
(positive downwards). The datum is `z = 0` for all models (reservoir and
wells), which is required for consistent energy exchange between them. Zero for
models that are not three-dimensional, consistent with how gravity is treated
in the flow equations.
"""
struct UnitPotentialEnergy <: ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::UnitPotentialEnergy, symb)
    # Note: The datum z = 0 is shared by all models (reservoir and wells), as
    # the cell centroids are given in the same global coordinates. This is
    # required for consistent exchange of energy between the models.
    cc = data_domain[:cell_centroids, Cells()]
    if size(cc, 1) == 3
        Φ = -gravity_constant.*vec(cc[3, :])
    else
        # Consistent with TwoPointGravityDifference: No gravity in 1D/2D
        Φ = zeros(size(cc, 2))
    end
    return Φ
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
    add_thermal_to_model!(model::MultiModel; energy_formulation = :thermal)
    add_thermal_to_model!(model; energy_formulation = :thermal)

Add energy conservation equation and thermal primary variable together with
standard set of parameters to existing flow model. Note that more complex models
require additional customization after this function call to get correct
results.

The `energy_formulation` keyword determines the conserved energy:
- `:thermal`: Conservation of thermal energy (internal energy of fluid and
  rock/well material), with advection of enthalpy and heat conduction.
- `:total`: Conservation of total energy, which in addition includes the
  gravitational potential energy of the fluid (see [`TotalEnergy`](@ref) and
  [`UnitPotentialEnergy`](@ref)). This accounts for the work done by gravity on
  the fluid, which can be significant in e.g. deep wells.

For a `MultiModel`, the same formulation is used for the reservoir and all
wells, which is required for consistent energy exchange between them.
"""
function add_thermal_to_model!(model::MultiModel; energy_formulation = :thermal)
    for (k, m) in pairs(model.models)
        if m.system isa MultiPhaseSystem
            add_thermal_to_model!(m, energy_formulation = energy_formulation)
        elseif m.system isa FacilitySystem
            add_thermal_to_facility!(m)
        end
    end
    return m
end

function add_thermal_to_model!(model; energy_formulation = :thermal)
    energy_formulation in (:thermal, :total) || throw(ArgumentError("energy_formulation must be :thermal or :total, was :$energy_formulation"))
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
    out = model.output_variables
    if energy_formulation == :total
        set_parameters!(model, UnitPotentialEnergy = UnitPotentialEnergy())
        set_secondary_variables!(model,
            PotentialEnergy = PotentialEnergy(),
            TotalEnergy = TotalEnergy()
        )
        model.equations[:energy_conservation] = ConservationLaw(disc, :TotalEnergy, 1)
        push!(out, :TotalEnergy)
    else
        model.equations[:energy_conservation] = ConservationLaw(disc, :TotalThermalEnergy, 1)
    end
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
    set_tabulated_internal_energy!(model, tab; regions = nothing)

Use tabulated fluid internal energy `tab(p, T)` (see
[`PressureTemperatureDependentInternalEnergy`](@ref)) for the reservoir and all
wells in `model`. The `ComponentHeatCapacity` variable is removed as it is no
longer used, so that it cannot be mixed with the tabulated values. Note that
`tab` must give the internal energy (not the enthalpy) and that all models must
use the same table so that the reference state is consistent.
"""
function set_tabulated_internal_energy!(model::MultiModel, tab; kwarg...)
    for (k, m) in pairs(model.models)
        if k == :Reservoir || model_or_domain_is_well(m)
            set_tabulated_internal_energy!(m, tab; kwarg...)
        end
    end
    return model
end

function set_tabulated_internal_energy!(model::SimulationModel, tab; regions = nothing)
    if !haskey(Jutul.get_secondary_variables(model), :FluidInternalEnergy)
        throw(ArgumentError("Model does not have FluidInternalEnergy, is it thermal?"))
    end
    if model_or_domain_is_well(model)
        # Wells are not partitioned by region
        regions = nothing
    end
    U = PressureTemperatureDependentInternalEnergy(tab, regions = regions)
    set_secondary_variables!(model, FluidInternalEnergy = U)
    Jutul.delete_variable!(model, :ComponentHeatCapacity)
    return model
end

"""
    model_energy_formulation(model)

Get the energy formulation (`:thermal` or `:total`) of a thermal model, or
`nothing` if the model has no energy equation.
"""
function model_energy_formulation(model::SimulationModel)
    eq = get(model.equations, :energy_conservation, nothing)
    if isnothing(eq)
        return nothing
    elseif eq isa ConservationLaw{:TotalEnergy}
        return :total
    else
        return :thermal
    end
end

function model_energy_formulation(model::MultiModel)
    return model_energy_formulation(reservoir_model(model))
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
