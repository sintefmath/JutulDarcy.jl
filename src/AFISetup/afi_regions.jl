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
        @assert !haskey(regions["family"], group) "Duplicate region family found: $group"
        gdict = Dict{String, Any}()
        regions["family"][group] = gdict

        for (i, k) in zip(indices, names)
            gdict[k] = i
        end
    end
    for prefix in ["Rock", "Fluid", "Equilibrium"]
        map_region!(regions, d, prefix)
    end
    # RockRegionMapping
    # FluidRegionMapping
    # EquilibriumRegionMapping
    return regions
end

function map_region!(regions, d::AFIInputFile, prefix)
    reg = Dict{String, Any}()
    regions[lowercase(prefix)] = reg
    recs = find_records(d, "$(prefix)RegionMapping", "IX", steps = true, model = true, once = false)
    for rec in recs
        @info "!" prefix rec.value
        tab = rec.value["table"]
        families = tab["RegionFamilyNames"]
        names = tab["RegionNames"]
        models = tab["ModelNames"]
        types = tab["MappingTypes"]
        for (i, fam, nm, mod, typ) in zip(1:length(families), families, names, models, types)
            if !haskey(reg, fam)
                reg[fam] = Dict{String, Any}()
            end
            fdict = reg[fam]
            if haskey(fdict, nm)
                @warn "Duplicate region mapping found for family $fam and name $nm. Overwriting previous entry."
            end
            fdict[nm] = (model = mod, type = String(typ))
        end
    end
end