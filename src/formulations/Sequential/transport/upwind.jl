function sort_tuple_indices(vals::NTuple{2, T}) where T
    a, b = vals
    if a < b
        ix = (1, 2)
    else
        ix = (2, 1)
    end
    return ix
end

function sort_tuple_indices(vals::NTuple{3, T}) where T
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

function phase_potential_upwind_fixed_flux(q, K, g::NTuple{N, T}, k_l::NTuple{N, V}, k_r::NTuple{N, V}) where {T, N, V}
    if N == 1
        flag = q < zero(q)
        flags = (flag,)
    else
        indices = sort_tuple_indices(g)
        if N == 2
            i1, i2 = indices
            Δg = g[i1] - g[i2]
            θ_1 = q + K*Δg*k_l[i2]
            θ_2 = q - K*Δg*k_r[i1]
            r = phase_potential_r_index(θ_1, θ_2)
        elseif N == 3
            i1, i2, i3 = indices
            Δg_12 = g[i1] - g[i2]
            Δg_13 = g[i1] - g[i3]
            Δg_23 = g[i2] - g[i3]
            # k_r = b
            θ_1 = q + K*(+Δg_12*k_l[i2] + Δg_13*k_l[i3])
            θ_2 = q + K*(-Δg_12*k_r[i1] + Δg_23*k_l[i3])
            θ_3 = q + K*(-Δg_13*k_r[i1] - Δg_23*k_r[i2])

            r = phase_potential_r_index(θ_1, θ_2, θ_3)
        else
            error("Not implemented for more than 3 phases")
        end
        # In paper: pick A (left) if l > r.
        # Our flags have opposite meaning: pick right if flag is true
        flags = indices .<= r
    end
    return flags
end

function phase_potential_upwind_potential_differences(V_t, T_f, G::NTuple{N, T}, left_mob, right_mob) where {N, T}
    flags = phase_potential_upwind_fixed_flux(V_t, T_f, G, left_mob, right_mob)

    function simple_upwind(l, r, flag::Bool)
        if flag
            v = r
        else
            v = l
        end
        return v
    end

    mob = map(simple_upwind, left_mob, right_mob, flags)
    mobT = sum(mob)

    F = 1.0/mobT
    if N == 2
        G_1, G_2 = G
        mob_1, mob_2 = mob

        ΔG_12 = G_1 - G_2
        dpot_1 = F*(V_t + T_f*ΔG_12*mob_2)
        dpot_2 = F*(V_t - T_f*ΔG_12*mob_1)
        phase_potential_differences = (dpot_1, dpot_2)
    else
        @assert N == 3
        G_1, G_2, G_3 = G
        mob_1, mob_2, mob_3 = mob

        ΔG_12 = G_1 - G_2
        ΔG_13 = G_1 - G_3
        ΔG_23 = G_2 - G_3
        dpot_1 = F*(V_t + T_f*ΔG_12*mob_2 + T_f*ΔG_13*mob_3)
        dpot_2 = F*(V_t - T_f*ΔG_12*mob_1 + T_f*ΔG_23*mob_3)
        dpot_3 = F*(V_t - T_f*ΔG_13*mob_1 - T_f*ΔG_23*mob_2)
        phase_potential_differences = (dpot_1, dpot_2, dpot_3)
    end
    return phase_potential_differences
end

function JutulDarcy.upwind(upw::SPU, F::AbstractArray, q::MultiPotential)
    new_pot = map(p -> JutulDarcy.upwind(upw, F, p), q.potentials)
    return MultiPotential(new_pot)
end

function JutulDarcy.phase_upwind(upw, m::AbstractMatrix, phase::Integer, q::MultiPotential)
    new_pot = map(p -> JutulDarcy.phase_upwind(upw, m, phase, p), q.potentials)
    return MultiPotential(new_pot)
end

function collapse_potentials(x::MultiPotential)
    return sum(x.potentials)
end

import Base.*

function (*)(q::MultiPotential, α::MultiPotential)
    new_pot = map((p, a) -> p*a, q.potentials, α.potentials)
    return MultiPotential(new_pot)
end

function (*)(q::MultiPotential, α::Number)
    new_pot = map(p -> p*α, q.potentials)
    return MultiPotential(new_pot)
end

function (*)(α::Number, q::MultiPotential)
    return q*α
end

import Base.+

function (+)(q::MultiPotential, p::MultiPotential)
    new_pot = map(+, q.potentials, p.potentials)
    return MultiPotential(new_pot)
end

function (+)(q::MultiPotential, α::Number)
    new_pot = map(p -> p + α, q.potentials)
    return MultiPotential(new_pot)
end

function (+)(α::Number, q::MultiPotential)
    return q + α
end

function Base.convert(t::Type{Float64}, q::MultiPotential)
    return convert(t, collapse_potentials(q))
end

using ForwardDiff
function Base.convert(t::Type{ForwardDiff.Dual{T, V, N}}, q::MultiPotential) where {T, V, N}
    return convert(t, collapse_potentials(q))
end
