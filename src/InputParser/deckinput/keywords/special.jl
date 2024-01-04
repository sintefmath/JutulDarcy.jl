function finish_current_section!(data, units, cfg, outer_data, ::Val{:EDIT})

end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:BOX})
    rec = read_record(f)
    tdims = [1];
    gdata = get_section(outer_data, :GRID)
    l, u = gdata["CURRENT_BOX"]
    il, jl, kl = l
    iu, ju, ku = u

    il, iu, jl, ju, kl, ku = parse_defaulted_line(rec, [il, iu, jl, ju, kl, ku])
    gdata["CURRENT_BOX"] = (lower = (il, jl, kl), upper = (iu, ju, ku))
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:ENDBOX})
    reset_current_box!(outer_data)
end

function parse_keyword!(data, outer_data, units, cfg, f, kval::Union{Val{:COPY}, Val{:ADD}, Val{:MULTIPLY}})
    k = unpack_val(kval)
    is_copy = k == :COPY
    rec = read_record(f)
    gdata = get_section(outer_data, :GRID)
    l, u = gdata["CURRENT_BOX"]
    dims = get_cartdims(outer_data)
    d = "Default"

    il = l[1]
    iu = u[1]
    jl = l[2]
    ju = u[2]
    kl = l[3]
    ku = u[3]
    if is_copy
        d_op = d
    else
        d_op = NaN
    end

    while length(rec) > 0
        parsed = parse_defaulted_line(rec, [d, d_op, il, iu, jl, ju, kl, ku])
        if is_copy
            dst = parsed[2]
            src = parsed[1]
            op = NaN
            @assert src != d "Source was defaulted for $k. rec = $rec"
        else
            dst = parsed[1]
            op = parsed[2]
            src = missing
            @assert op != d_op "Operator was defaulted for $k. rec = $rec"
        end
        @assert dst != d "Destination was defaulted for $k. rec = $rec"

        # Box can be kept.
        il = parsed[3]
        iu = parsed[4]
        jl = parsed[5]
        ju = parsed[6]
        kl = parsed[7]
        ku = parsed[8]
        if is_copy
            if !haskey(data, dst)
                T = eltype(data[src])
                data[dst] = zeros(T, dims)
            end
            apply_copy!(data[dst], data[src], get_box_indices(outer_data), dims)
        else
            if !haskey(data, dst)
                data[dst] = zeros(Float64, dims)
            end
            if k == :ADD
                # add is a const
                apply_add!(data[dst], op, get_box_indices(outer_data), dims)
            else
                # multiply is a const
                apply_multiply!(data[dst], op, get_box_indices(outer_data), dims)
            end
        end
        rec = read_record(f)
    end
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:OPERATE})
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

    parser_message(cfg, outer_data, "OPERATE", PARSER_MISSING_SUPPORT)

    while length(rec) > 0
        d = "Default"
        parsed = parse_defaulted_line(rec, [d, il, iu, jl, ju, kl, ku, d, d, NaN, NaN])
        src = parsed[1]
        op = parsed[8]
        @assert src != d "Source was defaulted? rec = $rec"
        @assert op != d "Operator was defaulted? rec = $rec"

        # Box can be kept.
        il = parsed[2]
        iu = parsed[3]
        jl = parsed[4]
        ju = parsed[5]
        kl = parsed[6]
        ku = parsed[7]

        op_prm1 = parsed[9]
        op_prm2 = parsed[10]
        op_prm3 = parsed[11]

        # @assert op_prm1 != "Default" "Operator parameter 1 was non-finite for OPERATE: $rec"
        # @assert isfinite(op_prm2) "Operator parameter 2 was non-finite for OPERATE: $rec"
        # @assert isfinite(op_prm3) "Operator parameter 3 was non-finite for OPERATE: $rec"
        # TODO: Implement operation
        # apply_copy!(data, dst, data[src], (il, iu), (jl, ju), (kl, ku), dims)
        rec = read_record(f)
    end
end

function apply_copy!(vals, src, IX, dims)
    I, J, K = IX
    vals[I, J, K] = src
end

function apply_add!(vals, src, IX, dims)
    I, J, K = IX
    vals[I, J, K] .+= src
end

function apply_multiply!(vals, src, IX, dims)
    I, J, K = IX
    vals[I, J, K] .*= src
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

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:EQUALS})
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
        constval = parsed[2]
        @assert dst != "Default"
        if haskey(data, dst)
            # Box can be kept.
            il = parsed[3]
            iu = parsed[4]
            jl = parsed[5]
            ju = parsed[6]
            kl = parsed[7]
            ku = parsed[8]
            apply_equals!(data[dst], constval, (il, iu), (jl, ju), (kl, ku), dims)
        else
            data[dst] = fill(constval, dims...)
        end
        rec = read_record(f)
    end
end

function apply_equals!(target, constval, I, J, K, dims)
    for i in I[1]:I[2]
        for j in J[1]:J[2]
            for k in K[1]:K[2]
                target[i, j, k] = constval
            end
        end
    end
end
