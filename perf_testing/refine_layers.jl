"""
Split layers of the two benchmark decks, retaining pillars and interpolating
ZCORN inside each original cell. This helper is deliberately scoped to these
decks, rather than being a general Eclipse deck refinement implementation.
"""
function refine_layers(input, factor::Int)
    factor >= 1 || throw(ArgumentError("refinement must be positive"))
    factor == 1 && return input
    data = deepcopy(input)
    grid = data["GRID"]
    dims = grid["cartDims"]
    nx, ny, nz = dims
    newdims = (nx, ny, nz * factor)
    @assert !any(haskey(grid, key) for key in ("NNC", "TRANX", "TRANY", "TRANZ", "AQUNUM", "AQUCON"))
    corners = reshape(grid["ZCORN"], 2nx, 2ny, 2nz)
    refined = similar(corners, 2nx, 2ny, 2nz * factor)
    for k in 1:nz, sub in 1:factor
        top = @view corners[:, :, 2k - 1]
        bottom = @view corners[:, :, 2k]
        child = (k - 1) * factor + sub
        @views refined[:, :, 2child - 1] .= top .+ ((sub - 1) / factor) .* (bottom .- top)
        @views refined[:, :, 2child] .= top .+ (sub / factor) .* (bottom .- top)
    end
    # Every cell field (including SATNUM, PVTNUM and endpoint scaling) is inherited.
    for section in values(data)
        section isa AbstractDict || continue
        for (key, value) in collect(pairs(section))
            value isa AbstractArray && size(value) == dims || continue
            section[key] = repeat(value; inner = (1, 1, factor))
            if key == "PORV"
                section[key] ./= factor
            elseif key in ("MULTZ", "MULTZ-")
                # A multiplier belongs to the original interface, not new faces.
                fill!(section[key], one(eltype(value)))
                offset = key == "MULTZ" ? factor : 1
                section[key][:, :, offset:factor:end] .= value
            end
        end
    end
    grid["ZCORN"] = vec(refined)
    grid["cartDims"] = newdims
    grid["CURRENT_BOX"] = (lower = (1, 1, 1), upper = newdims)
    data["RUNSPEC"]["DIMENS"] = collect(newdims)
    haskey(grid, "SPECGRID") && (grid["SPECGRID"][3] = newdims[3])
    haskey(grid, "MINPV") && (grid["MINPV"] /= factor)
    if haskey(grid, "FAULTS")
        for (name, faults) in grid["FAULTS"]
            refined_faults = empty(faults)
            for fault in faults
                if startswith(fault.direction, "Z") || startswith(fault.direction, "K")
                    for k in fault.k
                        child = endswith(fault.direction, "-") ? (k - 1) * factor + 1 : k * factor
                        push!(refined_faults, merge(fault, (k = child:child,)))
                    end
                else
                    k = ((first(fault.k) - 1) * factor + 1):(last(fault.k) * factor)
                    push!(refined_faults, merge(fault, (; k)))
                end
            end
            grid["FAULTS"][name] = refined_faults
        end
    end
    for step in data["SCHEDULE"]["STEPS"]
        @assert !any(haskey(step, key) for key in ("COMPSEGS", "COMPDATMD", "WPIMULT"))
        for completion in get(step, "COMPDAT", [])
            completion[4] = (completion[4] - 1) * factor + 1
            completion[5] *= factor
            # Preserve total supplied WI and Kh across the child completions.
            for index in (8, 10)
                isfinite(completion[index]) && completion[index] > 0 && (completion[index] /= factor)
            end
        end
        opens = get(step, "WELOPEN", nothing)
        isnothing(opens) && continue
        refined_opens = empty(opens)
        for record in opens
            @assert record[6] <= 0 && record[7] <= 0 "Connection-number WELOPEN is unsupported here"
            if record[5] > 0
                for sub in 1:factor
                    child = copy(record)
                    child[5] = (record[5] - 1) * factor + sub
                    push!(refined_opens, child)
                end
            else
                push!(refined_opens, record)
            end
        end
        step["WELOPEN"] = refined_opens
    end
    @assert count(grid["ACTNUM"]) == factor * count(input["GRID"]["ACTNUM"])
    return data
end
