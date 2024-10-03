function coarsen_reservoir_model(m::MultiModel, partition; functions = Dict(), well_arg = Dict(), kwarg...)
    reservoir = reservoir_domain(m)
    creservoir = coarsen_reservoir(reservoir, partition, functions = functions)

    wells = get_model_wells(m)
    cwells = []
    for (k, well) in pairs(wells)
        if haskey(well_arg, k)
            wkw = well_arg[k]
        else
            wkw = NamedTuple()
        end
        push!(cwells, coarsen_well(well, creservoir, reservoir, partition; wkw...))
    end
    split_wells = !haskey(m.models, :Facility)

    fine_reservoir_model = reservoir_model(m)
    sys = fine_reservoir_model.system

    coarse_model, coarse_parameters = setup_reservoir_model(creservoir, sys;
        wells = cwells,
        split_wells = split_wells,
        context = fine_reservoir_model.context,
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
    fine_keys = keys(m.models)
    @assert sort([coarse_keys...]) == sort([fine_keys...]) "Coarse keys $coarse_keys does not match fine keys $fine_keys"
    for (mkey, m) in pairs(coarse_model.models)
        if mkey == :Reservoir
            fmap = FiniteVolumeGlobalMap(subcells, Int[])
        elseif mkey in keys(wells)
            # A bit hackish to get this working for wells
            nnode = length(unique(m.domain.representation.perforations.self))
            fmap = FiniteVolumeGlobalMap(ones(Int, nnode), Int[])
        else
            fmap = TrivialGlobalMap()
        end
        for vartype in [:parameters, :primary, :secondary]
            vars = Jutul.get_variables_by_type(fine_reservoir_model, vartype)
            cvars = Jutul.get_variables_by_type(coarse_reservoir_model, vartype)
            for (k, var) in pairs(vars)
                cvars[k] = Jutul.subvariable(var, fmap)
            end
        end
    end
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
