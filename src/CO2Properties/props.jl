
function filter_table_to_temperature(sol, T)
    T_all, ix = get_column(sol, :T, index = true)
    Tu = unique(T_all)
    _, pos = findmin(abs.(Tu .- T))
    T_closest = Tu[pos]
    data = sol.data[T_all .== T_closest, :]
    return (data = data, header = copy(sol.header))
end

function co2brine_phase_property_for_T(data::Union{Tuple, AbstractVector}, name::Symbol, T::Float64)
    data = map(x -> filter_table_to_temperature(x, T), data)
    N = length(data)
    # Check that tables match up so that we can interpolate them together.
    d1 = first(data)
    for i in 2:N
        d = data[i]
        @assert size(d.data) == size(d1.data)
        @assert get_column(d, :p) == get_column(d1, :p)
        @assert get_column(d, :T) == get_column(d1, :T)
    end
    T_el = SVector{N, Float64}
    F = Vector{T_el}()
    p = get_column(d1, :p)
    sizehint!(F, length(p))
    tmp = zeros(N)
    for i in eachindex(p)
        for j in 1:N
            tmp[j] = get_column(data[j], name)[i]
        end
        push!(F, T_el(tmp))
    end
    return (p, F)
end

function co2brine_phase_property_table(data::NamedTuple, name, T = missing)
    return co2brine_phase_property_table((data, ), name, T)
end

function co2brine_phase_property_table(tables, name::Symbol, T = missing)
    choices = (:density, :viscosity, :H, :E, :cv, :cp, :phase_conductivity)
    name in choices || throw(ArgumentError("Property name must be one of $choices, was $name"))
    if ismissing(T)
        tab = first(tables)
        T = unique(get_column(tab, :T))
        p = unique(get_column(tab, :p))

        np = length(p)
        nt = length(T)
        property_val = zeros(SVector{length(tables), Float64}, np, nt)
        for (i, T_i) in enumerate(T)
            p_i, prop_i = co2brine_phase_property_for_T(tables, name, T_i)
            @assert p_i == p
            property_val[:, i] .= prop_i
        end
        @. T += 273.15
        property_eval = get_2d_interpolator(p, T, property_val)
    else
        p, property_val = co2brine_phase_property_for_T(tables, name, T - 273.15)
        property_eval = get_1d_interpolator(p, property_val, cap_endpoints = false)
    end
    return property_eval
end
