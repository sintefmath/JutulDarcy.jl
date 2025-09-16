function setup_region_map(d::AFIInputFile)
    region_defs = find_records(d, "RegionFamily", "IX", steps = true, model = true, once = false)
    regions = Dict{String, Any}(
        "family" => Dict{String, Any}()
    )

    for reg in region_defs
        reg = reg.value
        group = reg["group"]
        names = reg["RegionNames"]
        indices = get(reg, "RegionIndices", 1:length(names))
        if haskey(regions["family"], group)
            error("Duplicate region family found: $group")
        end
        gdict = Dict{String, Any}()
        regions["family"][group] = gdict

        for (i, k) in zip(indices, names)
            gdict[k] = i
        end
    end
    # RockRegionMapping
    # FluidRegionMapping
    # EquilibriumRegionMapping
    # ... others?
    for prefix in ["Rock", "Fluid", "Equilibrium"]
        map_region!(regions, d, prefix)
    end
    return regions
end

function map_region!(regions, d::AFIInputFile, prefix)
    reg = Dict{String, Any}()
    regions[lowercase(prefix)] = reg
    recs = find_records(d, "$(prefix)RegionMapping", "IX", steps = true, model = true, once = false)
    for rec in recs
        tab = rec.value["table"]
        families = tab["RegionFamilyNames"]
        names = tab["RegionNames"]
        models = tab["ModelNames"]
        types = tab["MappingTypes"]
        for (fam, nm, mod, typ) in zip(families, names, models, types)
            typ = String(typ)
            if !haskey(reg, typ)
                reg[typ] = Dict{String, Any}()
            end
            type_dict = reg[typ]
            if !haskey(type_dict, fam)
                type_dict[fam] = []
            end
            dest = type_dict[fam]
            push!(dest, (region = String(nm), model = String(mod)))
        end
    end
end
