module WellMerging

using Jutul
using DataStructures: OrderedDict

import ..JutulDarcy:
    SimpleWell, MultiSegmentWell, WellDomain, PerforationMask,
    add_thermal_to_model!, transfer_variables_and_parameters!,
    ReservoirFromWellFlowCT, ReservoirFromWellThermalCT,
    WellFromFacilityFlowCT, WellFromFacilityThermalCT,
    FacilityFromWellBottomHolePressureCT, FacilityFromSurfaceComponentRatesCT,
    FacilityFromWellTemperatureCT, FacilityFromWellEnthalpyCT,
    setup_reservoir_cross_terms!

"""
    merge_similar_wells(case::JutulCase)

Return a case where simple wells share one model and multisegment wells share
another. Each well remains independently controlled by its original name. The
merged domains carry the original names and node, face, and perforation ranges
in `multiwell`.
"""
function merge_similar_wells(case::JutulCase)
    original_model = case.model
    if !(original_model isa MultiModel)
        throw(ArgumentError("A reservoir MultiModel is required."))
    end

    simple_wells = unmerged_well_names(original_model, SimpleWell)
    multisegment_wells = unmerged_well_names(original_model, MultiSegmentWell)
    well_groups = (
        (:SimpleWells, simple_wells),
        (:MultiSegmentWells, multisegment_wells),
    )

    replacements = Dict{Symbol, Symbol}()
    merged_models = Dict{Symbol, Any}()
    for (merged_key, names) in well_groups
        if length(names) < 2
            continue
        end
        if haskey(original_model.models, merged_key)
            throw(ArgumentError("Model key $merged_key is already in use."))
        end

        merged_models[merged_key] = build_merged_well_model(
            original_model, names, merged_key)
        for name in names
            replacements[name] = merged_key
        end
    end

    if isempty(replacements)
        return case
    end

    new_models = OrderedDict{Symbol, Any}()
    source_indices = Dict{Symbol, Int}()
    for (index, (old_key, old_model)) in enumerate(pairs(original_model.models))
        new_key = get(replacements, old_key, old_key)
        if !haskey(new_models, new_key)
            if haskey(merged_models, new_key)
                new_models[new_key] = merged_models[new_key]
            else
                new_models[new_key] = old_model
            end
            source_indices[new_key] = index
        end
    end

    for pair in original_model.cross_terms
        if is_standard_well_cross_term(pair.cross_term)
            continue
        end
        target_is_merged = haskey(replacements, pair.target)
        source_is_merged = haskey(replacements, pair.source)
        if target_is_merged || source_is_merged
            throw(ArgumentError(
                "Cannot merge wells with an additional cross term attached to a well."))
        end
    end

    if isnothing(original_model.groups)
        new_groups = nothing
    else
        old_group_numbers = Int[]
        for key in keys(new_models)
            old_index = source_indices[key]
            push!(old_group_numbers, original_model.groups[old_index])
        end
        group_numbers = Dict{Int, Int}()
        new_groups = Int[]
        for old_group in old_group_numbers
            if !haskey(group_numbers, old_group)
                group_numbers[old_group] = length(group_numbers) + 1
            end
            push!(new_groups, group_numbers[old_group])
        end
    end

    modes = [original_model.group_execution[source_indices[key]]
        for key in keys(new_models)]
    merged_model = MultiModel(new_models;
        groups = new_groups, context = original_model.context,
        reduction = original_model.reduction,
        specialize_ad = original_model.specialize_ad,
        group_execution = modes)
    setup_reservoir_cross_terms!(merged_model)
    for pair in original_model.cross_terms
        if !is_standard_well_cross_term(pair.cross_term)
            push!(merged_model.cross_terms, pair)
        end
    end

    state0 = merge_case_storage(case.state0, merged_model, replacements)
    parameters = merge_case_storage(case.parameters, merged_model, replacements)
    if case.forces isa AbstractVector
        forces = [merge_case_forces(step_forces, merged_model, replacements)
            for step_forces in case.forces]
    else
        forces = merge_case_forces(case.forces, merged_model, replacements)
    end

    return JutulCase(merged_model, case.dt, forces, state0, parameters,
        case.input_data, case.termination_criterion, case.start_date)
end

function unmerged_well_names(model::MultiModel, well_type)
    names = Symbol[]
    for (name, submodel) in pairs(model.models)
        well = physical_representation(submodel)
        if well isa well_type && isnothing(well.multiwell)
            push!(names, name)
        end
    end
    return names
end

function build_merged_well_model(original_model, names, merged_key)
    well_models = [original_model.models[name] for name in names]
    template = first(well_models)
    system_type = typeof(template.system)
    same_systems = all(model -> typeof(model.system) == system_type, well_models)
    if !same_systems
        throw(ArgumentError("Wells being merged must have matching systems."))
    end

    domains = [model.data_domain for model in well_models]
    merged_domain = merge_well_domains(domains, merged_key)
    merged_model = SimulationModel(merged_domain, template.system;
        context = template.context, formulation = template.formulation,
        optimization_level = template.optimization_level)
    if haskey(template.primary_variables, :Temperature)
        add_thermal_to_model!(merged_model)
    end

    template_equations = Set(keys(template.equations))
    merged_equations = Set(keys(merged_model.equations))
    if merged_equations != template_equations
        throw(ArgumentError("Merging wells with additional equations is not supported."))
    end

    transfer_variables_and_parameters!(merged_model, template, check_type = false)
    empty!(merged_model.output_variables)
    append!(merged_model.output_variables, template.output_variables)
    return merged_model
end

is_standard_well_cross_term(ct) = ct isa Union{
    ReservoirFromWellFlowCT, ReservoirFromWellThermalCT,
    WellFromFacilityFlowCT, WellFromFacilityThermalCT,
    FacilityFromWellBottomHolePressureCT,
    FacilityFromSurfaceComponentRatesCT,
    FacilityFromWellTemperatureCT, FacilityFromWellEnthalpyCT}

function merge_well_domains(domains, name)
    wells = physical_representation.(domains)
    first_well = first(wells)
    names = Symbol[well.name for well in wells]
    cell_counts = number_of_cells.(wells)
    face_counts = number_of_faces.(wells)
    perforation_counts = [length(well.perforations.self) for well in wells]

    node_ranges = entity_ranges(cell_counts)
    face_ranges = entity_ranges(face_counts)
    perforation_ranges = entity_ranges(perforation_counts)
    top_nodes = [first(node_range) for node_range in node_ranges]
    info = (
        names = names,
        nodes = node_ranges,
        faces = face_ranges,
        perforations = perforation_ranges,
        top_nodes = top_nodes,
        domains = domains,
    )

    reservoir_cells = Int[]
    well_cells = Int[]
    for (well, node_range) in zip(wells, node_ranges)
        append!(reservoir_cells, well.perforations.reservoir)
        node_offset = first(node_range) - 1
        for node in well.perforations.self
            push!(well_cells, node + node_offset)
        end
    end
    perforations = (self = well_cells, reservoir = reservoir_cells)

    if first_well isa SimpleWell
        same_explicit_dp = all(well -> well.explicit_dp == first_well.explicit_dp, wells)
        if !same_explicit_dp
            throw(ArgumentError("Simple wells must have matching explicit_dp settings."))
        end
        merged_well = SimpleWell(perforations, first_well.surface, name,
            first_well.explicit_dp, info)
    else
        same_type = all(well -> well.type == first_well.type, wells)
        if !same_type
            throw(ArgumentError("Multisegment wells must have matching types."))
        end

        neighbor_blocks = []
        end_node_blocks = []
        segment_model_blocks = []
        for (well, node_range) in zip(wells, node_ranges)
            node_offset = first(node_range) - 1
            push!(neighbor_blocks, well.neighborship .+ node_offset)
            push!(end_node_blocks, well.end_nodes .+ node_offset)
            push!(segment_model_blocks, well.segment_models)
        end
        neighbors = reduce(hcat, neighbor_blocks)
        end_nodes = reduce(vcat, end_node_blocks)
        segment_models = reduce(vcat, segment_model_blocks)
        merged_well = MultiSegmentWell(first_well.type, sum(cell_counts),
            sum(face_counts), sum(perforation_counts), perforations, neighbors,
            end_nodes, first_well.surface, name, segment_models, info)
    end

    merged_domain = DataDomain(merged_well)
    data_keys = collect(keys(first(domains)))
    expected_keys = Set(data_keys)
    for domain in domains
        if Set(keys(domain)) != expected_keys
            throw(ArgumentError("Well data fields differ between wells being merged."))
        end
    end

    for key in data_keys
        entity = Jutul.associated_entity(first(domains), key)
        values = [domain[key] for domain in domains]
        if entity == NoEntity()
            same_values = all(isequal(first(values)), values)
            if !same_values
                throw(ArgumentError("Well property $key differs between wells."))
            end
            merged_domain[key, entity] = first(values)
        else
            merged_domain[key, entity] = cat(values...; dims = ndims(first(values)))
        end
    end
    return merged_domain
end

function entity_ranges(counts)
    ranges = UnitRange{Int}[]
    first_index = 1
    for count in counts
        last_index = first_index + count - 1
        push!(ranges, first_index:last_index)
        first_index = last_index + 1
    end
    return ranges
end

function merged_well_key(model::MultiModel, name::Symbol)
    if haskey(model.models, name)
        return name
    end
    for (key, submodel) in pairs(model.models)
        well = physical_representation(submodel)
        if !(well isa Union{SimpleWell, MultiSegmentWell})
            continue
        end
        if isnothing(well.multiwell)
            continue
        end
        if name in well.multiwell.names
            return key
        end
    end
    throw(KeyError(name))
end

function well_local_index(well::WellDomain, name::Symbol)
    if isnothing(well.multiwell)
        return 1
    end
    index = findfirst(isequal(name), well.multiwell.names)
    if isnothing(index)
        throw(KeyError(name))
    end
    return index
end

function well_top_node(well::WellDomain, name::Symbol)
    if isnothing(well.multiwell)
        return 1
    end
    index = well_local_index(well, name)
    return well.multiwell.top_nodes[index]
end

function well_perforations(well::WellDomain, name::Symbol)
    if isnothing(well.multiwell)
        return eachindex(well.perforations.self)
    end
    index = well_local_index(well, name)
    return well.multiwell.perforations[index]
end

function merge_entity_values(values)
    a = first(values)
    if a isa AbstractArray
        return cat(values...; dims = ndims(a))
    elseif a isa Number
        return collect(values)
    else
        throw(ArgumentError("Cannot merge well value of type $(typeof(a))."))
    end
end

function merge_case_storage(storage, model, replacements)
    if isnothing(storage)
        return nothing
    end

    merged_storage = OrderedDict{Symbol, Any}()
    for key in keys(model.models)
        if key in values(replacements)
            well = physical_representation(model.models[key])
            entries = [storage[name] for name in well.multiwell.names]
            merged_entry = OrderedDict{Symbol, Any}()
            fields = collect(keys(first(entries)))
            expected_fields = Set(fields)
            same_fields = all(entry -> Set(keys(entry)) == expected_fields, entries)
            if !same_fields
                throw(ArgumentError("Well state or parameter fields differ between wells."))
            end
            for field in fields
                field_values = [entry[field] for entry in entries]
                merged_entry[field] = merge_entity_values(field_values)
            end
            merged_storage[key] = merged_entry
        else
            merged_storage[key] = storage[key]
        end
    end
    return merged_storage
end

function merge_case_forces(forces, model, replacements)
    merged_forces = OrderedDict{Symbol, Any}()
    for key in keys(model.models)
        if key in values(replacements)
            well = physical_representation(model.models[key])
            masks = [forces[name].mask for name in well.multiwell.names]
            if all(isnothing, masks)
                merged_forces[key] = (mask = nothing,)
            else
                mask_values = []
                for (mask, perforation_range) in zip(masks, well.multiwell.perforations)
                    if isnothing(mask)
                        push!(mask_values, ones(length(perforation_range)))
                    else
                        push!(mask_values, mask.values)
                    end
                end
                combined_values = reduce(vcat, mask_values)
                merged_forces[key] = (mask = PerforationMask(combined_values),)
            end
        else
            merged_forces[key] = forces[key]
        end
    end
    return merged_forces
end

end # module WellMerging
