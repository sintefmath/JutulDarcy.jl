"""
    coarsen_reservoir_model(fine_model::MultiModel, partition; functions = Dict(), well_arg = Dict(), kwarg...)

Coarsens a fine reservoir model based on the provided partitioning scheme.

# Arguments
- `fine_model::MultiModel`: The fine-scale reservoir model to be coarsened.
- `partition`: The vector used to coarsen the model.
- `functions`: A dictionary of functions to be applied during the coarsening process. Defaults to an empty dictionary.
- `well_arg`: A dictionary of well arguments to be considered during coarsening. Defaults to an empty dictionary.
- `kwarg...`: Additional keyword arguments that are passed onto (`setup_reservoir_model`)[@ref].

# Returns
- A coarsened reservoir model.
"""
function coarsen_reservoir_model(fine_model::MultiModel, partition; functions = Dict(), well_arg = Dict(), kwarg...)
    reservoir = reservoir_domain(fine_model)
    creservoir = coarsen_reservoir(reservoir, partition, functions = functions)

    wells = get_model_wells(fine_model)
    cwells = []
    for (k, well) in pairs(wells)
        if haskey(well_arg, k)
            wkw = well_arg[k]
        else
            wkw = NamedTuple()
        end
        push!(cwells, coarsen_well(well, creservoir, reservoir, partition; wkw...))
    end
    split_wells = !haskey(fine_model.models, :Facility)

    fine_reservoir_model = reservoir_model(fine_model)
    sys = fine_reservoir_model.system

    thermal = haskey(Jutul.get_primary_variables(fine_reservoir_model), :Temperature)

    coarse_model, = setup_reservoir_model(creservoir, sys;
        wells = cwells,
        split_wells = split_wells,
        thermal = thermal,
        kwarg...
    )
    coarse_reservoir_model = reservoir_model(coarse_model)
    # Variables etc.
    ncoarse = maximum(partition)
    subcells = zeros(Int, ncoarse)
    for i in 1:ncoarse
        subcells[i] = findfirst(isequal(i), partition)
    end
    coarse_keys = keys(coarse_model.models)
    fine_keys = keys(fine_model.models)
    @assert sort([coarse_keys...]) == sort([fine_keys...]) "Coarse keys $coarse_keys does not match fine keys $fine_keys"
    for mkey in keys(coarse_model.models)
        cm = coarse_model.models[mkey]
        fm = fine_model.models[mkey]
        if mkey == :Reservoir
            fmap = FiniteVolumeGlobalMap(subcells, Int[])
            fmap::FiniteVolumeGlobalMap
        elseif mkey in keys(wells)
            # A bit hackish to get this working for wells
            nnode = length(unique(cm.domain.representation.perforations.self))
            fmap = FiniteVolumeGlobalMap(ones(Int, nnode), Int[])
        else
            fmap = TrivialGlobalMap()
        end
        for vartype in [:parameters, :primary, :secondary]
            vars = Jutul.get_variables_by_type(fm, vartype)
            cvars = Jutul.get_variables_by_type(cm, vartype)
            empty!(cvars)
            for (k, var) in pairs(vars)
                cvars[k] = Jutul.subvariable(var, fmap)
            end
        end
    end
    coarse_parameters = setup_parameters(coarse_model)
    return (coarse_model, coarse_parameters)
end

function coarsen_reservoir(D::DataDomain, partition; functions = Dict())
    if !haskey(functions, :permeability)
        functions[:permeability] = Jutul.CoarsenByHarmonicAverage()
    end
    if !haskey(functions, :porosity)
        functions[:porosity] = Jutul.CoarsenByVolumeAverage()
    end
    return Jutul.coarsen_data_domain(D, partition, functions = functions)
end

function coarsen_well(well, creservoir::DataDomain, reservoir::DataDomain, partition; kwarg...)
    cells = unique(partition[well.perforations.reservoir])
    is_simple = well isa SimpleWell
    if is_simple
        d = well.reference_depth
    else
        d = well.top.reference_depth
    end
    return setup_well(creservoir, cells;
        name = well.name,
        reference_depth = d,
        simple_well = is_simple,
        kwarg...
    )
end

struct CoarsenByPoreVolume <: Jutul.AbstractCoarseningFunction
    fine_pv::Vector{Float64}
end

function Jutul.inner_apply_coarsening_function(finevals, fine_indices, op::CoarsenByPoreVolume, coarse, fine, row, name, entity)
    subvols = op.fine_pv[fine_indices]
    return sum(finevals.*subvols)/sum(subvols)
end

"""
    coarsen_reservoir_state(coarse_model, fine_model, fine_state0; functions = Dict(), default = missing)

Coarsens the reservoir state from a fine model to a coarse model.

# Arguments
- `coarse_model`: The coarse model to which the state will be coarsened.
- `fine_model`: The fine model from which the state will be coarsened.
- `fine_state0`: The initial state of the fine model.
- `functions`: A dictionary of functions to apply during the coarsening process. Defaults to an empty dictionary.
- `default`: The default value to use if no function is provided for a specific operation. Defaults to `missing`.

# Returns
- The coarsened state of the reservoir.
"""
function coarsen_reservoir_state(coarse_model, fine_model, fine_state0; functions = Dict(), default = missing)
    coarse_state0 = Dict{Symbol, Any}()

    if haskey(fine_state0, :Reservoir)
        # We ignore the wells, better to re-initialize
        fine_state0 = fine_state0[:Reservoir]
    end
    coarse_rmodel = reservoir_model(coarse_model)
    fine_rmodel = reservoir_model(fine_model)

    coarse_reservoir = reservoir_domain(coarse_rmodel)
    fine_reservoir = reservoir_domain(fine_rmodel)

    if ismissing(default)
        pv = pore_volume(fine_reservoir)
        default = CoarsenByPoreVolume(pv)
    end

    ncoarse = number_of_cells(coarse_reservoir)
    for (k, v) in pairs(fine_state0)
        if associated_entity(fine_rmodel[k]) == Cells() && eltype(v)<:AbstractFloat
            if v isa AbstractVector
                coarseval = zeros(ncoarse)
            else
                coarseval = zeros(size(v, 1), ncoarse)
            end
            f = get(functions, k, default)
            coarse_state0[k] = Jutul.apply_coarsening_function!(coarseval, v, f, coarse_reservoir, fine_reservoir, k, Cells())
        end
    end
    return setup_reservoir_state(coarse_model, coarse_state0)
end

function partition_reservoir(case::JutulCase, coarsedim, method = missing; kwarg...)
    return partition_reservoir(case.model, coarsedim, method; parameters = case.parameters, kwarg...)
end

"""
    partition_reservoir(model::JutulModel, coarsedim::Union{Tuple, Int}, method = missing)

Partition the reservoir model into coarser grids.

# Arguments
- `model::JutulModel`: The reservoir model to be partitioned.
- `coarsedim::Union{Tuple, Int}`: The dimensions for the coarser grid. Can be a
  tuple specifying dimensions or an integer to specify the total number of
  desired coarse blocks.
- `method`: Optional. The method to use for partitioning. Defaults to `missing`.
- `partitioner_conn_type`: Optional. The type of connection to use for the
  partition. Can be :trans, :logtrans or :unit.

# Returns
- A partitioned version of the reservoir model.

"""
function partition_reservoir(model::JutulModel, coarsedim::Union{Tuple, Int}, method = missing;
        parameters = missing,
        wells_in_single_block = false,
        partitioner_conn_type = :trans
    )
    domain = model |> reservoir_model |> reservoir_domain
    mesh = physical_representation(domain)

    if coarsedim isa Int
        method = :metis
        coarsedim > 0 || throw(ArgumentError("Number of blocks must be positive, was $coarsedim."))
    else
        length(coarsedim) == dim(mesh) || throw(ArgumentError("coarsedim argument must be tuple of equal length to dimension of grid (e.g. (5, 5, 5) for 3D)"))
        if ismissing(method)
            method = :centroids
        else
            method in (:ijk, :centroids) || throw(ArgumentError("When coarsedim is a tuple, method must be either :ijk or :centroids"))
        end
    end

    if method in (:ijk, :centroids)
        p = Jutul.cartesian_partition(domain, coarsedim, method)
    else
        if method == :metis
            partitioner = Jutul.MetisPartitioner()
        elseif method == :linear
            partitioner = Jutul.LinearPartitioner()
        elseif method == :kahypar
            partitioner = Jutul.KaHyParPartitioner()
        else
            partitioner = method
        end
        partitioner::Jutul.JutulPartitioner
        if ismissing(parameters)
            parameters = setup_parameters(model)
        end
        N, T, well_groups = partitioner_input(model, parameters, conn = partitioner_conn_type)
        if wells_in_single_block
            groups = well_groups
        else
            groups = Vector{Vector{Int}}()
        end
        p = Jutul.partition_hypergraph(N, coarsedim, partitioner, groups = groups)
        p = Int64.(p)
    end
    p = Jutul.process_partition(mesh, p)
    p = Jutul.compress_partition(p)
    return p
end

export coarsen_reservoir_case


"""
    coarsen_reservoir_case(case, coarsedim; kwargs...)

Coarsens the given reservoir case to the specified dimensions.

# Arguments
- `case`: The reservoir case to be coarsened.
- `coarsedim`: The target dimensions for the coarsened reservoir.

# Keyword Arguments
- `method`: The method to use for partitioning. Defaults to `missing`.
- `partitioner_arg`: A named tuple of arguments to be passed to the partitioner.
- `setup_arg`: A named tuple of arguments to be passed to `coarsen_reservoir_model`.
- `state_arg`: A named tuple of arguments to be passed to the state coarsening function.

# Returns
- A coarsened version of the reservoir case.

"""
function coarsen_reservoir_case(case, coarsedim;
        method = missing,
        partitioner_arg = NamedTuple(),
        setup_arg = NamedTuple(),
        state_arg = NamedTuple()
    )
    if coarsedim isa Vector
        p = coarsedim
    else
        p = partition_reservoir(case, coarsedim, method; partitioner_arg...)
    end
    (; model, forces, dt, parameters, state0) = case
    coarse_model, coarse_parameters = coarsen_reservoir_model(model, p; setup_arg...)
    coarse_state0 = coarsen_reservoir_state(coarse_model, model, state0; state_arg...)
    coarse_forces = deepcopy(forces)
    coarse_dt = deepcopy(dt)
    return JutulCase(coarse_model, coarse_dt, coarse_forces, parameters = coarse_parameters, state0 = coarse_state0)
end
