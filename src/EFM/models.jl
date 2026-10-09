"""
    setup_embedded_fracture_model(matrix, fractures, system; representation=:coupled, subtract_fracture_volume=false, kwarg...)

Set up an EFM reservoir model from a preprocessed fracture domain.

* `representation=:coupled`: separate `:Reservoir` and `:Fractures` models with
  the same conservative matrix–fracture cross terms as DFM. Supports immiscible,
  single-phase, and `:geothermal` systems.
* `representation=:single`: one reservoir model with fracture cells appended
  after all matrix cells. Supports systems accepted by `setup_reservoir_model`.

Interior matrix connections remain open. Face-coincident fractures replace the
covered fraction of direct matrix connections. Optional volume subtraction
accounts for fracture storage without modifying either input domain. Wells
keep their original matrix perforations; explicit fracture well completions
can be assigned to appended cell indices in the single representation.
"""
function setup_embedded_fracture_model(matrix::DataDomain, fractures::DataDomain, system;
        representation = :coupled,
        subtract_fracture_volume = false,
        extra_out = false,
        wells = [],
        kwarg...
    )
    mesh = physical_representation(fractures)
    mesh isa EmbeddedFractureMesh || throw(ArgumentError("fractures must be created with embedded_fracture_domain"))
    physical_representation(matrix) === mesh.matrix_mesh || throw(ArgumentError("EFM domains must share the original matrix mesh"))
    if representation == :single
        domain = combined_fracture_domain(matrix, fractures; subtract_fracture_volume = subtract_fracture_volume)
        model = setup_reservoir_model(domain, system; wells = wells, kwarg...)
        if JutulDarcy.model_is_thermal(reservoir_model(model))
            set_parameters!(reservoir_model(model), FluidThermalConductivities = EFMFluidThermalConductivities())
        end
        if extra_out
            return model, setup_parameters(model)
        end
        return model
    elseif representation != :coupled
        throw(ArgumentError("representation must be :coupled or :single"))
    end
    system isa Union{ImmiscibleSystem, SinglePhaseSystem} || system == :geothermal || throw(ArgumentError("coupled EFM supports immiscible, single-phase, and :geothermal systems; use representation=:single for other systems"))
    matrix = embedded_matrix_domain(matrix, fractures; subtract_fracture_volume = subtract_fracture_volume)
    model = setup_reservoir_model(matrix, system; wells = wells, kwarg...)
    fracture_model = reservoir_model(setup_reservoir_model(fractures, system; kwarg...))
    fmc = FractureMatrixConnection()
    fracture_model.domain.entities[fmc] = count_entities(fractures, fmc)
    if JutulDarcy.model_is_thermal(fracture_model)
        set_parameters!(fracture_model, FluidThermalConductivities = EFMFluidThermalConductivities())
    end
    model = JutulDarcy.add_fracture_model(model, fracture_model)
    JutulDarcy.add_fracture_cross_terms!(model, fractures)
    if extra_out
        return model, setup_parameters(model)
    end
    return model
end

function JutulDarcy.setup_fractured_reservoir_model(matrix::DataDomain, fractures::DataDomain{<:EmbeddedFractureMesh}, system::Union{JutulSystem, Symbol}; kwarg...)
    return setup_embedded_fracture_model(matrix, fractures, system; kwarg...)
end
struct EFMFluidThermalConductivities <: VectorVariables end
Jutul.associated_entity(::EFMFluidThermalConductivities) = Faces()
Jutul.minimum_value(::EFMFluidThermalConductivities) = 0.0
Jutul.variable_scale(::EFMFluidThermalConductivities) = 1.0
Jutul.values_per_entity(model, ::EFMFluidThermalConductivities) = number_of_phases(model.system)

function Jutul.default_parameter_values(domain, model, ::EFMFluidThermalConductivities, name)
    T = domain[:fluid_thermal_conductivities]
    nph = number_of_phases(model.system)
    if T isa AbstractVector
        return repeat(T', nph, 1)
    elseif size(T, 1) == 1
        return repeat(T, nph, 1)
    end
    size(T, 1) == nph || throw(ArgumentError("fluid thermal conductivities must have one row per phase"))
    return copy(T)
end
