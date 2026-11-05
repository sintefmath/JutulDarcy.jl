# Translate deck groups to the flat, fixed membership supported by WellGroup.
# Production and each injection phase get separate allocation groups, even when
# they share a WELSPECS group name.
function deck_group_control(rec, injection = false)
    group = rec[1]
    mode = rec[ifelse(injection, 3, 2)]
    if injection
        mode in ("NONE", "FLD", "RATE") || throw(ArgumentError("GCONINJE $group: $mode group control is not supported (only surface RATE)."))
        if mode == "NONE"
            any(isfinite, rec[4:7]) && throw(ArgumentError("GCONINJE $group: constraints with NONE control are not supported."))
            return nothing
        end
        rec[9] == 0.0 && isempty(rec[10]) || throw(ArgumentError("GCONINJE $group: explicit group guide rates are not supported."))
        if mode == "FLD"
            any(isfinite, rec[4:7]) && throw(ArgumentError("GCONINJE $group: subordinate group constraints are not supported."))
            return nothing
        end
        phase = rec[2]
        if phase == "WATER" || phase == "WAT"
            Target = GroupSurfaceWaterRateTarget
        elseif phase == "OIL"
            Target = GroupSurfaceOilRateTarget
        elseif phase == "GAS"
            Target = GroupSurfaceGasRateTarget
        else
            throw(ArgumentError("GCONINJE $group: unsupported injection phase $phase."))
        end
        any(isfinite, rec[5:7]) && throw(ArgumentError("GCONINJE $group: reservoir rate, reinjection and voidage constraints are not supported."))
        isfinite(rec[4]) && rec[4] >= 0 || throw(ArgumentError("GCONINJE $group RATE requires a nonnegative surface target."))
        target = Target(rec[4])
        limits = as_limit(target)
    else
        modes = ("ORAT", "WRAT", "GRAT", "LRAT")
        types = (GroupSurfaceOilRateTarget, GroupSurfaceWaterRateTarget,
            GroupSurfaceGasRateTarget, GroupSurfaceLiquidRateTarget)
        mode in (modes..., "NONE", "FLD") || throw(ArgumentError("GCONPROD $group: $mode group control is not supported."))
        if mode == "NONE"
            (any(isfinite, rec[3:6]) || isfinite(rec[14])) &&
                throw(ArgumentError("GCONPROD $group: constraints with NONE control are not supported."))
            return nothing
        end
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
            # The parser uses NONE for both omitted and explicit procedures.
            # Those mean different things when the general procedure is RATE.
            # Require an explicit procedure for additional phase constraints
            # until that distinction is retained by the parser.
            if i != ix && i != 1 && isfinite(rec[i + 2]) && rec[7] == "RATE" && rec[9 + i] == "NONE"
                throw(ArgumentError("GCONPROD $group: the parser cannot distinguish defaulted and explicit NONE exceed procedures. Specify RATE explicitly for additional phase constraints."))
            end
            if i == 1 || rec[9 + i] == "NONE"
                procedure = rec[7]
            else
                procedure = rec[9 + i]
            end
            if procedure == "RATE" && isfinite(rec[i + 2])
                rec[i + 2] >= 0 || throw(ArgumentError("GCONPROD $group has a negative rate constraint."))
                limits = merge(limits, as_limit(Target(-rec[i + 2])))
            end
        end
    end
    return (; target, limits)
end

function parse_deck_group_steps(schedule)
    steps = schedule["STEPS"]
    groups = Dict{Symbol, Vector{Symbol}}()
    if !any(step -> haskey(step, "GCONPROD") || haskey(step, "GCONINJE"), steps)
        # Guide-rate records and group trees do not affect ordinary well
        # controls when the deck requests no group constraints.
        return (; steps = [Dict{String, Any}() for _ in steps], groups)
    end
    parents = Dict{String, String}()
    memberships = Dict(name => spec.group for (name, spec) in schedule["WELSPECS"])
    production = Dict{String, Any}()
    injection = Dict{Tuple{String, String}, Any}()
    kinds = Dict{String, String}()
    participation = Dict{String, Bool}()
    owners = Dict{String, Symbol}()
    all_specs = []
    has_guide_formula = false
    for step in steps
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
            phase = ifelse(rec[2] == "WAT", "WATER", rec[2])
            injection[(rec[1], phase)] = rec
        end
        for rec in get(step, "WGRUPCON", [])
            rec[2] in ("YES", "NO") || throw(ArgumentError("Invalid WGRUPCON participation flag $(rec[2])."))
            if rec[2] == "YES"
                rec[3] <= 0 && rec[5] == 1.0 || throw(ArgumentError("WGRUPCON: explicit well guide rates are not supported."))
            end
            if haskey(participation, rec[1]) && participation[rec[1]] != (rec[2] == "YES")
                throw(ArgumentError("WGRUPCON: changing group participation during the schedule is not supported."))
            end
            participation[rec[1]] = rec[2] == "YES"
        end
        for key in ("WCONPROD", "WCONHIST", "WCONINJE", "WCONINJ", "WCONINJH")
            for rec in get(step, key, [])
                kind = ifelse(key in ("WCONPROD", "WCONHIST"), "PRODUCTION", rec[2])
                kinds[rec[1]] = ifelse(kind == "WAT", "WATER", kind)
            end
        end
        specs = Dict{String, Any}()
        step_owners = Dict{String, Symbol}()
        for (well, kind) in kinds
            is_injection = kind != "PRODUCTION"
            group = memberships[well]
            seen = Set{String}()
            selected = nothing
            while true
                group in seen && throw(ArgumentError("Cycle in GRUPTREE at $group."))
                push!(seen, group)
                if is_injection
                    rec = get(injection, (group, kind), nothing)
                else
                    rec = get(production, group, nothing)
                end
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
                step_owners[well] = owner
                members = get!(groups, owner, Symbol[])
                Symbol(well) in members || push!(members, Symbol(well))
                # Opted-out wells keep their individual controls and limits,
                # but their rates still consume part of the group target.
                if get(participation, well, true)
                    specs[well] = selected
                end
            end
        end
        for (well, owner) in owners
            if haskey(kinds, well) && !haskey(step_owners, well) && any(isequal(owner), values(step_owners))
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
