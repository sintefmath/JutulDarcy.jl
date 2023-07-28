function unpack_val(::Val{X}) where X
    return X
end

const DECK_SPLIT_REGEX = r"[ \t,]+"

function read_record(f; fix = true)
    split_lines = Vector{String}()
    active = true
    while !eof(f) && active
        line = strip(readline(f))
        if !startswith(line, "--")
            if endswith(line, '/')
                line = strip(rstrip(line, '/'))
                active = false
            end
            if length(line) > 0
                push!(split_lines, line)
            end
        end
    end
    return split_lines
end

function keyword_start(line)
    if isnothing(match(r"^\s*--", line))
        m = match(r"\w+", line)
        if m === nothing
            return nothing
        else
            return Symbol(uppercase(m.match))
        end
    else
        return nothing
    end
end

function parse_defaulted_group(f, defaults)
    out = []
    line = read_record(f)
    while length(line) > 0
        parsed = parse_defaulted_line(line, defaults)
        push!(out, parsed)
        line = read_record(f)
    end
    return out
end
##
function parse_defaulted_line(lines::String, defaults)
    return parse_defaulted_line([lines], defaults)
end

function parse_defaulted_line(lines, defaults)
    out = similar(defaults, 0)
    sizehint!(out, length(defaults))
    pos = 1
    for line in lines
        lsplit = split(strip(line), DECK_SPLIT_REGEX)
        for s in lsplit
            if length(s) == 0
                continue
            end
            if occursin('*', s) && !startswith(s, '\'') # Could be inside a string for wildcard matching
                if s == "*"
                    num_defaulted = 1
                else
                    num_defaulted = Parsers.parse(Int, match(r"\d+\*", s).match[1:end-1])
                end
                for i in 1:num_defaulted
                    push!(out, defaults[pos])
                    pos += 1
                end
            else
                default = defaults[pos]
                if default isa String
                    converted = strip(s, [' ', '\''])
                else
                    T = typeof(default)
                    converted = Parsers.parse(T, s)
                end
                push!(out, converted)
                pos += 1
            end
        end
    end
    n = length(defaults)
    if pos < n
        for i in pos:n
            push!(out, defaults[i])
        end
    end
    return out
end

##

function parse_deck_matrix(f, T = Float64)
    # TODO: This is probably a bad way to do large numerical datasets.
    rec = read_record(f)
    split_lines = preprocess_delim_records(rec)
    data = Vector{T}()
    n = -1
    for seg in split_lines
        m = length(seg)
        if m == 0
            continue
        elseif n == -1
            n = m
        else
            @assert m == n "Expected $n was $m"
        end
        for d in seg
            push!(data, parse(T, d))
        end
    end
    return reshape(data, n, length(data) ÷ n)'
end

function preprocess_delim_records(split_lines)
    # Strip end whitespace
    split_lines = map(strip, split_lines)
    # Remove comments
    filter!(x -> !startswith(x, "--"), split_lines)
    # Split into entries (could be whitespace + a comma anywhere in between)
    split_rec = map(x -> split(x, r"\s*,?\s+"), split_lines)
    # Remove entries
    for recs in split_rec
        filter!(x -> length(x) > 0, recs)
    end
    return split_rec
end

function parse_deck_vector(f, T = Float64)
    # TODO: Speed up.
    rec = read_record(f)
    record_lines = preprocess_delim_records(rec)
    n = length(record_lines)
    out = Vector{T}()
    sizehint!(out, n)
    for split_rec in record_lines
        for el in split_rec
            if occursin('*', el)
                n_rep, el = split(el, '*')
                n_rep = parse(Int, n_rep)
            else
                n_rep = 1
            end
            for i in 1:n_rep
                parsed = parse(T, el)::T
                push!(out, parsed)
            end
        end
    end
    return out
end

function skip_record(f)
    rec = read_record(f)
    while length(rec) > 0
        rec = read_record(f)
    end
end

function skip_records(f, n)
    for i = 1:n
        rec = read_record(f)
    end
end

function parse_grid_vector(f, dims, T = Float64)
    v = parse_deck_vector(f, T)
    return reshape(v, dims)
end

function parse_saturation_table(f, outer_data)
    ns = number_of_tables(outer_data, :saturation)
    return parse_region_matrix_table(f, ns)
end

function parse_dead_pvt_table(f, outer_data)
    np = number_of_tables(outer_data, :pvt)
    return parse_region_matrix_table(f, np)
end

function parse_live_pvt_table(f, outer_data)
    nreg = number_of_tables(outer_data, :pvt)
    out = []
    for i = 1:nreg
        current = Vector{Vector{Float64}}()
        while true
            next = parse_deck_vector(f)
            if length(next) == 0
                break
            end
            push!(current, next)
        end
        push!(out, restructure_pvt_table(current))
    end
    return out
end

function restructure_pvt_table(tab)
    nvals_per_rec = 3
    function record_length(x)
        # Actual number of records: 1 key value + nrec*N entries. Return N.
        return (length(x) - 1) ÷ nvals_per_rec
    end
    @assert record_length(last(tab)) > 1
    nrecords = length(tab)
    keys = map(first, tab)
    current = 1
    for tab_ix in eachindex(tab)
        rec = tab[tab_ix]
        interpolate_missing_usat!(tab, tab_ix, record_length, nvals_per_rec)
    end
    # Generate final table
    ntab = sum(record_length, tab)
    data = zeros(ntab, nvals_per_rec)
    for tab_ix in eachindex(tab)
        rec = tab[tab_ix]
        for i in 1:record_length(rec)
            for j in 1:nvals_per_rec
                linear_ix = (i-1)*nvals_per_rec + j + 1
                data[current, j] = rec[linear_ix]
            end
            current += 1
        end
    end

    # Generate pos
    pos = Int[1]
    sizehint!(pos, nrecords+1)
    for rec in tab
        push!(pos, pos[end] + record_length(rec))
    end
    return Dict("data" => data, "key" => keys, "pos" => pos)
end

function interpolate_missing_usat!(tab, tab_ix, record_length, nvals_per_rec)
    rec = tab[tab_ix]
    if record_length(rec) == 1
        @assert nvals_per_rec == 3
        next_rec = missing
        for j in (tab_ix):length(tab)
            if record_length(tab[j]) > 1
                next_rec = tab[j]
                break
            end
        end
        @assert record_length(rec) == 1
        next_rec_length = record_length(next_rec)
        sizehint!(rec, 1 + nvals_per_rec*next_rec_length)

        get_index(major, minor) = nvals_per_rec*(major-1) + minor + 1
        pressure(x, idx) = x[get_index(idx, 1)]
        B(x, idx) = x[get_index(idx, 2)]
        viscosity(x, idx) = x[get_index(idx, 3)]

        function constant_comp_interp(F, F_l, F_r)
            # So that dF/dp * F = constant over the pair of points extrapolated from F
            w = 2.0*(F_l - F_r)/(F_l + F_r)
            return F*(1.0 + w/2.0)/(1.0 - w/2.0)
        end
        @assert !ismissing(next_rec) "Final table must be saturated."

        for idx in 2:next_rec_length
            # Each of these gets added as new unsaturated points
            p_0 = pressure(rec, idx - 1)
            p_l = pressure(next_rec, idx - 1)
            p_r = pressure(next_rec, idx)

            mu_0 = viscosity(rec, idx - 1)
            mu_l = viscosity(next_rec, idx - 1)
            mu_r = viscosity(next_rec, idx)

            B_0 = B(rec, idx - 1)
            B_l = B(next_rec, idx - 1)
            B_r = B(next_rec, idx)

            p_next = p_0 + p_r - p_l
            B_next = constant_comp_interp(B_0, B_l, B_r)
            mu_next = constant_comp_interp(mu_0, mu_l, mu_r)

            push!(rec, p_next)
            push!(rec, B_next)
            push!(rec, mu_next)
        end
    end
    return tab
end

function parse_region_matrix_table(f, nreg)
    out = []
    for i = 1:nreg
        push!(out, parse_deck_matrix(f))
    end
    return out
end

function parse_keyword!(data, outer_data, units, f, ::Val{T}) where T
    # Do nothing
    error("Unhandled keyword $T encountered.")
end

function parse_keyword!(data, outer_data, units, f, ::Val{:EQUIL})
    n = number_of_tables(outer_data, :equil)
    def = [0.0, NaN, 0.0, 0.0, 0.0, 0.0, 0, 0, 0]
    out = []
    for i = 1:n
        rec = read_record(f)
        result = parse_defaulted_line(rec, def)
        push!(out, result)
    end
    return out
end


function next_keyword!(f)
    m = nothing
    while isnothing(m) && !eof(f)
        line = readline(f)
        m = keyword_start(line)
    end
    return m
end

function number_of_tables(outer_data, t::Symbol)
    td = outer_data["RUNSPEC"]["TABDIMS"]
    if t == :saturation
        return td[1]
    elseif t == :pvt
        return td[2]
    elseif t == :equil
        return outer_data["RUNSPEC"]["EQLDIMS"][1]
    else
        error(":$t is not known")
    end
end

function clean_include_path(basedir, include_file_name)
    include_file_name = strip(include_file_name)
    include_file_name = replace(include_file_name, "./" => "")
    include_file_name = replace(include_file_name, "'" => "")
    include_path = joinpath(basedir, include_file_name)
    return include_path
end

function get_section(outer_data, name::Symbol)
    s = "$name"
    is_sched = name == :SCHEDULE
    T = OrderedDict{String, Any}
    if is_sched
        if !haskey(outer_data, s)
            outer_data[s] = [T()]
        end
        out = outer_data[s][end]
    else
        if !haskey(outer_data, s)
            outer_data[s] = T()
        end
        out = outer_data[s]
    end
    return out
end

function new_section(outer_data, name::Symbol)
    data = get_section(outer_data, name)
    return data
end
