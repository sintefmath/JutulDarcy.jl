
function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:COPY})
    rec = read_record(f)
    gdata = get_section(outer_data, :GRID)
    l, u = gdata["CURRENT_BOX"]
    dims = get_cartdims(outer_data)

    il = l[1]
    iu = u[1]
    jl = l[2]
    ju = u[2]
    kl = l[3]
    ku = u[3]

    while length(rec) > 0
        d = "Default"
        parsed = parse_defaulted_line(rec, [d, d, il, iu, jl, ju, kl, ku])
        src = parsed[1]
        dst = parsed[2]
        @assert src != "Default" "Source was defaulted? rec = $rec"
        @assert dst != "Default" "Destination was defaulted? rec = $rec"

        # Box can be kept.
        il = parsed[3]
        iu = parsed[4]
        jl = parsed[5]
        ju = parsed[6]
        kl = parsed[7]
        ku = parsed[8]
        apply_copy!(data, dst, data[src], (il, iu), (jl, ju), (kl, ku), dims)
        rec = read_record(f)
    end
end

function apply_copy!(data, dst, src, I, J, K, dims)
    if haskey(data, dst)
        error("Not implemented")
    else
        @assert I[1] == J[1] == K[1] == 1
        @assert I[2] == dims[1]
        @assert J[2] == dims[2]
        @assert K[2] == dims[3]

        data[dst] = copy(src)
    end
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MULTIPLY})
    # TODO: Merge shared code with COPY
    rec = read_record(f)
    l, u = outer_data["GRID"]["CURRENT_BOX"]
    dims = outer_data["GRID"]["cartDims"]

    il = l[1]
    iu = u[1]
    jl = l[2]
    ju = u[2]
    kl = l[3]
    ku = u[3]

    while length(rec) > 0
        d = "Default"
        parsed = parse_defaulted_line(rec, [d, 1.0, il, iu, jl, ju, kl, ku])
        dst = parsed[1]
        factor = parsed[2]
        @assert dst != "Default"

        # Box can be kept.
        il = parsed[3]
        iu = parsed[4]
        jl = parsed[5]
        ju = parsed[6]
        kl = parsed[7]
        ku = parsed[8]
        @assert haskey(data, dst) "Unable to apply MULTIPLY to non-declared field $dst"
        apply_multiply!(data[dst], factor, (il, iu), (jl, ju), (kl, ku), dims)
        rec = read_record(f)
    end
end

function apply_multiply!(target::AbstractVector, factor, I, J, K, dims)
    apply_multiply!(reshape(target, dims), factor, I, J, K, dims)
end


function apply_multiply!(target, factor, I, J, K, dims)
    for i in I[1]:I[2]
        for j in J[1]:J[2]
            for k in K[1]:K[2]
                target[i, j, k] *= factor
            end
        end
    end
end
