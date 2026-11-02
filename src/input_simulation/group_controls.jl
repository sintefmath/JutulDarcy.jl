# Translate deck groups to the flat, fixed membership supported by WellGroup.
# Production and each injection phase get separate allocation groups, even when
# they share a WELSPECS group name.
function deck_group_control(rec, injection = false)
    group = rec[1]
    mode = rec[injection ? 3 : 2]
    if injection
        mode in ("NONE", "FLD", "RATE") || throw(ArgumentError("GCONINJE $group: $mode group control is not supported (only surface RATE)."))
        mode == "NONE" && return nothing
        rec[9] == 0.0 && isempty(rec[10]) || throw(ArgumentError("GCONINJE $group: explicit group guide rates are not supported."))
        if mode == "FLD"
            any(isfinite, rec[4:7]) && throw(ArgumentError("GCONINJE $group: subordinate group constraints are not supported."))
            return nothing
        end
        phase = rec[2]
        Target = phase == "WATER" ? GroupSurfaceWaterRateTarget :
            phase == "OIL" ? GroupSurfaceOilRateTarget :
            phase == "GAS" ? GroupSurfaceGasRateTarget : nothing
        isnothing(Target) && throw(ArgumentError("GCONINJE $group: unsupported injection phase $phase."))
        any(isfinite, rec[5:7]) && throw(ArgumentError("GCONINJE $group: reservoir rate, reinjection and voidage constraints are not supported."))
        isfinite(rec[4]) && rec[4] >= 0 || throw(ArgumentError("GCONINJE $group RATE requires a nonnegative surface target."))
        target = Target(rec[4])
        limits = as_limit(target)
    else
        modes = ("ORAT", "WRAT", "GRAT", "LRAT")
        types = (GroupSurfaceOilRateTarget, GroupSurfaceWaterRateTarget,
            GroupSurfaceGasRateTarget, GroupSurfaceLiquidRateTarget)
        mode in (modes..., "NONE", "FLD") || throw(ArgumentError("GCONPROD $group: $mode group control is not supported."))
        mode == "NONE" && return nothing
        !isfinite(rec[9]) && isempty(rec[10]) || throw(ArgumentError("GCONPROD $group: explicit group guide rates are not supported."))
        all(p -> p in ("NONE", "RATE"), rec[[7, 11, 12, 13]]) ||
            throw(ArgumentError("GCONPROD $group: only RATE exceed procedures are supported."))
        isfinite(rec[14]) && throw(ArgumentError("GCONPROD $group: reservoir-volume group constraints are not supported."))
        if mode == "FLD"
            any(isfinite, rec[3:6]) && throw(ArgumentError("GCONPROD $group: subordinate group constraints are not supported."))
            return nothing
        end
        ix = findfirst(isequal(mode), modes)
        rate = rec[ix + 2]
        isfinite(rate) && rate >= 0 || throw(ArgumentError("GCONPROD $group $mode requires a nonnegative surface target."))
        target = types[ix](-rate)
        limits = as_limit(target)
        for (i, Target) in enumerate(types)
            procedure = i == 1 ? rec[7] : rec[9 + i] == "NONE" ? rec[7] : rec[9 + i]
            if procedure == "RATE" && isfinite(rec[i + 2])
                rec[i + 2] >= 0 || throw(ArgumentError("GCONPROD $group has a negative rate constraint."))
                limits = merge(limits, as_limit(Target(-rec[i + 2])))
            end
        end
    end
    return (; target, limits)
end

function parse_deck_group_steps(schedule)
    parents = Dict{String, String}()
    memberships = Dict(name => spec.group for (name, spec) in schedule["WELSPECS"])
    production = Dict{String, Any}()
    injection = Dict{Tuple{String, String}, Any}()
    kinds = Dict{String, String}()
    participation = Dict{String, Bool}()
    groups = Dict{Symbol, Vector{Symbol}}()
    owners = Dict{String, Symbol}()
    all_specs = []
    has_guide_formula = false
    for step in schedule["STEPS"]
        has_guide_formula |= haskey(step, "GUIDERAT")
        for rec in get(step, "GRUPTREE", [])
            parents[rec[1]] = rec[2]
        end
        for rec in get(step, "WELSPECS", [])
            memberships[rec[1]] = rec[2]
        end
        for rec in get(step, "GCONPROD", [])
            occursin(r"[*?]", rec[1]) && throw(ArgumentError("GCONPROD: wildcard group names are not supported."))
            production[rec[1]] = rec
        end
        for rec in get(step, "GCONINJE", [])
            occursin(r"[*?]", rec[1]) && throw(ArgumentError("GCONINJE: wildcard group names are not supported."))
            injection[(rec[1], rec[2])] = rec
        end
        for rec in get(step, "WGRUPCON", [])
            rec[2] in ("YES", "NO") || throw(ArgumentError("Invalid WGRUPCON participation flag $(rec[2])."))
            rec[3] < 0 && rec[5] == 1.0 || throw(ArgumentError("WGRUPCON: explicit well guide rates are not supported."))
            if haskey(participation, rec[1]) && participation[rec[1]] != (rec[2] == "YES")
                throw(ArgumentError("WGRUPCON: changing group participation during the schedule is not supported."))
            end
            participation[rec[1]] = rec[2] == "YES"
        end
        for key in ("WCONPROD", "WCONHIST", "WCONINJE", "WCONINJ", "WCONINJH")
            for rec in get(step, key, [])
                kinds[rec[1]] = key in ("WCONPROD", "WCONHIST") ? "PRODUCTION" : rec[2]
            end
        end
        specs = Dict{String, Any}()
        for (well, kind) in kinds
            get!(participation, well, true) || continue
            is_injection = kind != "PRODUCTION"
            group = memberships[well]
            seen = Set{String}()
            selected = nothing
            while true
                group in seen && throw(ArgumentError("Cycle in GRUPTREE at $group."))
                push!(seen, group)
                rec = is_injection ? get(injection, (group, kind), nothing) : get(production, group, nothing)
                if !isnothing(rec)
                    control = deck_group_control(rec, is_injection)
                    if !isnothing(control)
                        isnothing(selected) || throw(ArgumentError("Well $well has constraints at multiple group levels; nested group constraints are not supported."))
                        selected = (; control..., group = Symbol(group, "_", kind))
                    end
                    rec[8] == "NO" && break
                end
                group == "FIELD" && break
                group = get(parents, group, "FIELD")
            end
            if !isnothing(selected)
                owner = selected.group
                if haskey(owners, well) && owners[well] != owner
                    throw(ArgumentError("Well $well changes allocation group during the schedule; dynamic group membership is not supported."))
                end
                owners[well] = owner
                members = get!(groups, owner, Symbol[])
                Symbol(well) in members || push!(members, Symbol(well))
                specs[well] = selected
            end
        end
        for (well, owner) in owners
            if haskey(kinds, well) && !haskey(specs, well) && any(spec -> spec.group == owner, values(specs))
                throw(ArgumentError("Well $well leaves an active allocation group; dynamic group membership is not supported."))
            end
        end
        if !isempty(specs) && has_guide_formula
            throw(ArgumentError("GUIDERAT formulas are not supported by deck group control."))
        end
        push!(all_specs, specs)
    end
    return (; steps = all_specs, groups)
end

function apply_deck_group_limits(limits, spec)
    isnothing(spec) && return limits
    isnothing(limits) && return spec.limits
    return merge(limits, spec.limits)
end
