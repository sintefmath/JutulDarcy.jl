function K_values_for_T(sol, T)
    sol = filter_table_to_temperature(sol, T)
    p = get_column(sol, :p)
    K = SVector{2, Float64}[]

    myclamp(x) = clamp(x, 1e-20, 1.0)
    for i in eachindex(p)
        y_h2o = myclamp(get_column(sol, :y_h2o)[i])
        x_co2 = myclamp(get_column(sol, :x_co2)[i])

        x_h2o = myclamp(1.0 - x_co2)
        y_co2 = myclamp(1.0 - y_h2o)
        # Kx = y
        K_h2o = y_h2o/x_h2o
        K_co2 = y_co2/x_co2
        if K_h2o < 0 || !isfinite(K_h2o)
            @warn "Bad K_value:" K_h2o y_h2o x_h2o y_co2 x_co2
        end
        if K_co2 < 0 || !isfinite(K_co2)
            @warn "Bad K_value:" K_co2 y_h2o x_h2o y_co2 x_co2
        end
        K_i = SVector{2, Float64}(K_h2o, K_co2)
        push!(K, K_i)
    end
    return (p, K)
end

function co2brine_K_values(sol, T = missing)
    if ismissing(T)
        T = unique(get_column(sol, :T))
        p = unique(get_column(sol, :p))

        np = length(p)
        nt = length(T)

        K = zeros(SVector{2, Float64}, np, nt)
        for (i, T_i) in enumerate(T)
            p_i, K_i = K_values_for_T(sol, T_i)
            @assert p_i == p
            K[:, i] .= K_i
        end
        # To Kelvin
        @. T += 273.15
        K_eval = get_2d_interpolator(p, T, K)
        dep = :pT
        out = KValueWrapper(K_eval, dependence = :pT)
    else
        p, K = K_values_for_T(sol, T - 273.153)
        K_eval = get_1d_interpolator(p, K, cap_endpoints = true)
        dep = :p
    end
    return KValueWrapper(K_eval, dependence = dep)
end

