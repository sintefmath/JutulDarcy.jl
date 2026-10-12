"""
    setup_embedded_fracture_model(matrix, fractures, system; unified=false, extend_wells=false, subtract_fracture_volume=false, kwarg...)

Set up an EFM reservoir model from a preprocessed fracture domain.

* `unified=false`: separate `:Reservoir` and `:Fractures` models with
  the same conservative matrix–fracture cross terms as DFM. Supports immiscible,
  single-phase, and `:geothermal` systems.
* `unified=true`: one reservoir model with fracture cells appended
  after all matrix cells. Supports systems accepted by `setup_reservoir_model`.
  Its domain has an integer `:is_fracture` cell field: `0` for matrix cells and
  `1` for fracture cells.

Interior matrix connections remain open. Face-coincident fractures replace the
covered fraction of direct matrix connections. Optional volume subtraction
accounts for fracture storage without modifying either input domain. Wells
keep their original matrix perforations; explicit fracture well completions
can be assigned to appended cell indices in the unified model.
Set `extend_wells=true` to add fracture completions where multisegment well
segments intersect clipped fracture cells, in either formulation.
"""
function setup_embedded_fracture_model(matrix::DataDomain, fractures::DataDomain, system;
        unified::Bool = false,
        subtract_fracture_volume = false,
        extra_out = false,
        wells = [],
        extend_wells = false,
        kwarg...
    )
    mesh = physical_representation(fractures)
    mesh isa EmbeddedFractureMesh || throw(ArgumentError("fractures must be created with embedded_fracture_domain"))
    physical_representation(matrix) === mesh.matrix_mesh || throw(ArgumentError("EFM domains must share the original matrix mesh"))
    if extend_wells
        wells = extend_embedded_wells(wells, fractures)
    end
    if unified
        domain = combined_fracture_domain(matrix, fractures; subtract_fracture_volume = subtract_fracture_volume)
        wells = [unified_fracture_well(well, number_of_cells(matrix)) for well in wells]
        model = setup_reservoir_model(domain, system; wells = wells, kwarg...)
        if JutulDarcy.model_is_thermal(reservoir_model(model))
            set_parameters!(reservoir_model(model), FluidThermalConductivities = EFMFluidThermalConductivities())
        end
    else
        system isa Union{ImmiscibleSystem, SinglePhaseSystem} || system == :geothermal || throw(ArgumentError("coupled EFM supports immiscible, single-phase, and :geothermal systems; use unified=true for other systems"))
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
    end
    out = model
    if extra_out
        out = (model, setup_parameters(model))
    end
    return out
end

function JutulDarcy.setup_fractured_reservoir_model(matrix::DataDomain, fractures::DataDomain{<:EmbeddedFractureMesh}, system::Union{JutulSystem, Symbol}; kwarg...)
    return setup_embedded_fracture_model(matrix, fractures, system; kwarg...)
end

const EmbeddedFractureObjects = Union{AbstractEmbeddedFracture, AbstractVector, Tuple}

function setup_embedded_fracture_model(matrix::DataDomain, objects::EmbeddedFractureObjects, system;
        fracture_arg = NamedTuple(), kwarg...)
    fractures = embedded_fracture_domain(matrix, objects; fracture_arg...)
    return setup_embedded_fracture_model(matrix, fractures, system; kwarg...)
end

"""
    setup_fractured_reservoir_model(matrix, objects, system; method=:dfm, unified=false, fracture_arg=NamedTuple(), wells=[], kwarg...)

Build either cut-cell DFM (`method=:dfm`, default) or EFM (`method=:efm`) from
a vector of fracture objects. Any vector element type is accepted; entries are
validated and grouped in a `FractureCollection`. For either method, `unified=true` appends
fracture cells to the reservoir model; `unified=false` uses separate matrix and
fracture models. Preprocessing options and property
overrides are passed in `fracture_arg`. DFM well perforations are distributed over
their cut-cell descendants, preserving each completion's total well index.
`extend_wells=false` keeps matrix completions; set it to `true` to add
fracture completions to multisegment wells for either method.
Both methods default to `subtract_fracture_volume=true` for matching storage.
The original domain and well data are unchanged. Use `discrete_fracture_domain`
when the original-cell mapping is needed for transferring states or forces.
"""
function JutulDarcy.setup_fractured_reservoir_model(matrix::DataDomain, objects::EmbeddedFractureObjects,
        system::Union{JutulSystem, Symbol};
        method = :dfm, unified::Bool = false,
        fracture_arg = NamedTuple(), wells = [], extend_wells = false,
        subtract_fracture_volume = true, extra_out = false, kwarg...)
    method in (:efm, :dfm) || throw(ArgumentError("method must be :efm or :dfm"))
    objects = FractureCollection(objects)
    if method == :efm
        model = setup_embedded_fracture_model(matrix, objects, system;
            unified = unified, fracture_arg = fracture_arg,
            wells = wells, extend_wells = extend_wells,
            subtract_fracture_volume = subtract_fracture_volume, kwarg...)
    else
        domains = discrete_fracture_domain(matrix, objects; fracture_arg...)
        remapped = [cut_fracture_well(well, domains.matrix, domains.info[:cell_index]) for well in wells]
        model = setup_fractured_reservoir_model(domains.matrix, domains.fractures, system;
            unified = unified, wells = remapped, extend_wells = extend_wells,
            subtract_fracture_volume = subtract_fracture_volume, kwarg...)
    end
    out = model
    if extra_out
        out = (model, setup_parameters(model))
    end
    return out
end

function cut_fracture_well(well::DataDomain, matrix, parents)
    original = physical_representation(well)
    original isa Union{JutulDarcy.SimpleWell, JutulDarcy.MultiSegmentWell} || throw(ArgumentError("DFM well transfer requires a SimpleWell or MultiSegmentWell"))
    isnothing(original.multiwell) || throw(ArgumentError("transfer individual wells before combining them into a multiwell"))
    children = [Int[] for _ in 1:maximum(parents)]
    for (cell, parent) in enumerate(parents)
        push!(children[parent], cell)
    end
    cells, self, indices = Int[], Int[], Int[]
    weights = eltype(matrix[:volumes])[]
    for (perf, parent) in enumerate(original.perforations.reservoir)
        descendants = children[parent]
        volume = sum(matrix[:volumes][descendants])
        for cell in descendants
            push!(cells, cell)
            push!(self, original.perforations.self[perf])
            push!(indices, perf)
            push!(weights, matrix[:volumes][cell]/volume)
        end
    end
    perforations = (self = self, reservoir = cells)
    w = deepcopy(original)
    if w isa JutulDarcy.SimpleWell
        representation = JutulDarcy.SimpleWell(perforations, w.surface, w.name, w.explicit_dp, w.multiwell)
    else
        representation = JutulDarcy.MultiSegmentWell(w.type, w.num_nodes, w.num_segments,
            length(cells), perforations, w.neighborship, w.end_nodes, w.surface,
            w.name, w.segment_models, w.include_inertia, w.multiwell)
    end
    transferred = DataDomain(representation)
    for (key, (value, entity)) in pairs(well)
        if entity == JutulDarcy.Perforations()
            mapped = fracture_cell_subset(value, indices)
            if key == :well_index
                values = Jutul.default_parameter_values(well, nothing, JutulDarcy.WellIndices(), :WellIndices)
                mapped = values[indices].*weights
            elseif key == :thermal_well_index && haskey(well, :thermal_conductivity)
                values = Jutul.default_parameter_values(well, nothing, JutulDarcy.WellIndicesThermal(), :WellIndicesThermal)
                mapped = values[indices].*weights
            elseif key == :perforation_centroids
                mapped = matrix[:cell_centroids][:, cells]
            elseif key == :cell_dims
                mesh = physical_representation(matrix)
                mapped = [JutulDarcy.peaceman_cell_dims(mesh, cell) for cell in cells]
            end
            transferred[key, entity] = mapped
        else
            transferred[key, entity] = deepcopy(value)
        end
    end
    return transferred
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
