"""
    kozeny_carman_permeability(porosity, constant = 1.0si_unit(:darcy); bounds = (1e-20, Inf))
    kozeny_carman_permeability(porosity, sphericity, grain_diameter, shape_factor, tortuosity; bounds = (1e-20, Inf))

Compute permeability from porosity using the Kozeny-Carman relation

``k = C \\frac{ϕ^3}{(1 - ϕ)^2}``

clamped to `bounds`. The constant `C` can either be given directly, or computed
from the grain properties as ``C = ψ^2 d^2 / (36 k_0 τ)``, where ``ψ`` is the
sphericity, ``d`` the grain diameter, ``k_0`` the shape factor and ``τ`` the
tortuosity. `porosity` can be a scalar or an array.
"""
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

function kozeny_carman_permeability(porosity, sphericity, grain_diameter, shape_factor, tortuosity; kwargs...)

    ψ = sphericity
    d = grain_diameter
    k0 = shape_factor
    τ = tortuosity
    constant = (ψ^2*d^2)/(36*k0*τ)

    return kozeny_carman_permeability(porosity, constant; kwargs...)

end
