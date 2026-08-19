function kozeny_carman_permeability(
    porosity, constant::Real = 1.0si_unit(:darcy);
    bounds::Tuple{<:Real, <:Real} = (1e-20, Inf)
    )

    constant > 0 || throw(ArgumentError("constant must be positive."))
    kmin, kmax = bounds
    0.0 < kmin <= kmax || throw(ArgumentError("bounds must satisfy 0 < min <= max."))

    ϕ = porosity
    C = constant
    k = C * (ϕ.^3) ./ ((1 .- ϕ).^2)

    return clamp.(k, kmin, kmax)

end

function kozeny_carman_permeability(porosity, sphericity, grain_diameter, shape_factor, tourtosity; kwargs...)

    ψ = sphericity
    d = grain_diameter
    k0 = shape_factor
    τ = tourtosity
    constant = (ψ^2*d^2)/(36*k0*τ)
    
    return kozeny_carman_permeability(porosity, constant; kwargs...)

end
