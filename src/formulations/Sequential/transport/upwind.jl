function sort_tuple(vals::NTuple{2, T}) where T
    a, b = vals
    if a < b
        ix = (1, 2)
    else
        ix = (2, 1)
    end
    return ix
end

function sort_tuple(vals::NTuple{3, T}) where T
    v1, v2, v3 = vals
    if v1 < v2
        if v1 < v3
            if v2 < v3
                ix = (1, 2, 3)
            else
                ix = (1, 3, 2)
            end
        else
            ix = (3, 1, 2)
        end 
    else
        if v2 < v3
            if v1 < v3
                ix = (2, 1, 3)
            else
                ix = (2, 3, 1)
            end
        else
            ix = (3, 2, 1)
        end
    end
    return ix
end

function phase_potential_r_index(thetas...)
    r = 0
    for i in length(thetas):-1:1
        theta = thetas[i]
        if theta <= 0
            r = i
            break
        end
    end
    return r
end

function phase_potential_upwind_fixed_flux(q, K, g::NTuple{N, T}, k_l::NTuple{N, V}, k_r::NTuple{N, V}, debug::Bool = false) where {T, N, V}
    debug = false
    # debug = true
    if N == 1
        flags = q < zero(q)
    else
        vals = q .+ K.*g
        indices = sort_tuple(vals)
        if N == 2
            i1, i2 = indices
            Δg = g[i1] - g[i2]
            θ_1 = q + K*Δg*k_r[i2]
            θ_2 = q - K*Δg*k_l[i1]
            pos_1 = θ_1 >= 0
            pos_2 = θ_2 >= 0
            r = phase_potential_r_index(θ_1, θ_2)
            if debug
                @info "" q θ_1 θ_2 r k_l k_r g indices
            end
        elseif N == 3
            i1, i2, i3 = indices
            Δg_12 = g[i1] - g[i2]
            Δg_13 = g[i1] - g[i3]
            Δg_23 = g[i2] - g[i3]
            θ_1 = q + K*(Δg_12*k_r[i2] + Δg_13*k_r[i3])
            θ_2 = q + K*(-Δg_12*k_l[i1] + Δg_23*k_r[i3])
            θ_3 = q + K*(-Δg_13*k_l[i1] - Δg_23*k_l[i2])

            pos_1 = θ_1 >= 0
            pos_2 = θ_2 >= 0
            pos_3 = θ_3 >= 0

            r = phase_potential_r_index(θ_1, θ_2, θ_3)
            if debug
                @info "" q θ_1 θ_2 θ_3 r k_l k_r g indices
            end
        else
            error("Not implemented for more than 3 phases")
        end
        flags = indices .<= r
    end
    return flags
end
