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
        IJK = get_box_indices(outer_data, il, iu, jl, ju, kl, ku)
        if is_copy
            if !haskey(data, dst)
                T = eltype(data[src])
                data[dst] = zeros(T, dims)
            end
            apply_copy!(data[dst], data[src], IJK, dims)
        else
            if !haskey(data, dst)
                data[dst] = zeros(Float64, dims)
            end
            if k == :ADD
                # add is a const
                apply_add!(data[dst], op, IJK, dims)
            else
                # multiply is a const
                apply_multiply!(data[dst], op, IJK, dims)
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

    while length(rec) > 0
        d = "Default"
        parsed = parse_defaulted_line(rec, [d, il, iu, jl, ju, kl, ku, d, d, NaN, NaN])
        target = parsed[1]
        op = parsed[8]
        source = parsed[9]
        @assert target != d "Target was defaulted? rec = $rec"
        @assert op != d "Operator was defaulted? rec = $rec"
        @assert source != d "Source was defaulted? rec = $rec"

        # Box can be kept.
        il = parsed[2]
        iu = parsed[3]
        jl = parsed[4]
        ju = parsed[5]
        kl = parsed[6]
        ku = parsed[7]

        op_prm1 = parsed[10]
        op_prm2 = parsed[11]

        IJK = get_box_indices(outer_data, il, iu, jl, ju, kl, ku)
        operation_target = get_operation_section(outer_data, target)
        operation_source = get_operation_section(outer_data, source)
        if ismissing(operation_target)
            data[target] = zeros(dims)
            operation_target = data[target]
        end
        apply_operate!(operation_target, operation_source, IJK, op, op_prm1, op_prm2)
        # On to the next one.
        rec = read_record(f)
    end
end

function apply_operate!(target, source, IJK, op, prm1, prm2)
    I, J, K = IJK
    op = lowercase(op)
    if op == "multx"
        F = (t, s) -> prm1*s
    elseif op == "addx"
        F = (t, s) -> s + prm1
    elseif op == "multa"
        F = (t, s) -> prm1*s + prm2
    elseif op == "abs"
        F = (t, s) -> abs(s)
    elseif op ==  "loge"
        F = (t, s) -> ln(s)
    elseif op ==  "log10"
        F = (t, s) -> log10(s)
    elseif op ==  "slog"
        F = (t, s) -> 10^(prm1 + prm2*s)
    elseif op ==  "poly"
        F = (t, s) -> t + prm1*s^prm2
    elseif op ==  "inv"
        F = (t, s) -> 1.0/s
    elseif op == "multiply"
        F = (t, s) -> t*s
    elseif op == "multp"
        F = (t, s) -> prm1*s^prm2
    elseif op == "minlim"
        F = (t, s) -> max(prm1, s)
    elseif op == "maxlim"
        F = (t, s) -> min(prm1, s)
    elseif op == "copy"
        F = (t, s) -> s
    else
        error("OPERATE option $(uppercase(op)) is not implemented.")
    end
    for i in I
        for j in J
            for k in K
                target[i, j, k] = F(target[i, j, k], source[i, j, k])
            end
        end
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

function get_operation_section(outer_data, kw)
    for (k, data) in pairs(outer_data)
        if data isa AbstractDict && haskey(data, kw)
            return data[kw]
        end
    end
    return missing
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MULTIPLY})
    # TODO: Merge shared code with COPY
    rec = read_record(f)
    grid = outer_data["GRID"]
    l, u = grid["CURRENT_BOX"]
    dims = grid["cartDims"]

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
        target = get_operation_section(outer_data, dst)
        if ismissing(target)
            if dst == "PORV" && haskey(grid, "PORO")
                # TODO: Bit of a hack
                target = grid["PORO"]
            else
                throw(ArgumentError("Unable to apply MULTIPLY to non-declared field $dst"))
            end
        end
        apply_multiply!(target, factor, (il, iu), (jl, ju), (kl, ku), dims)
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
        target = get_operation_section(outer_data, dst)
        if ismissing(target)
            # TODO: Different keywords go in different spots...
            data[dst] = fill(constval, dims...)
        else
            # Box can be kept.
            il = parsed[3]
            iu = parsed[4]
            jl = parsed[5]
            ju = parsed[6]
            kl = parsed[7]
            ku = parsed[8]
            apply_equals!(target, constval, (il, iu), (jl, ju), (kl, ku), dims)
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

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MAXVALUE})
    edit_apply_clamping!(f, outer_data, min)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MINVALUE})
    edit_apply_clamping!(f, outer_data, max)
end

function edit_apply_clamping!(f, outer_data, FUNCTION)
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
        lim = parsed[2]
        @assert dst != "Default"
        target = get_operation_section(outer_data, dst)
        # Box can be kept.
        il = parsed[3]
        iu = parsed[4]
        jl = parsed[5]
        ju = parsed[6]
        kl = parsed[7]
        ku = parsed[8]
        apply_equals!(target, lim, (il, iu), (jl, ju), (kl, ku), dims)
        rec = read_record(f)
    end
end

function apply_minmax!(target, lim, I, J, K, dims, F)
    for i in I[1]:I[2]
        for j in J[1]:J[2]
            for k in K[1]:K[2]
                target[i, j, k] = F(target[i, j, k], lim)
            end
        end
    end
end
