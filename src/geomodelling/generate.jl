function _standard_normal_generator(dims...)
    return randn(dims...)
end

"""
    generate_perm_poro(dims, box_lengths = (1.0, 1.0, 1.0); kwargs...)

Generate porosity and permeability fields from a latent Gaussian draw.

The real implementation lives in the GeoStats extension when that package is
loaded. This stub keeps the public API available in the main package.
"""
function generate_perm_poro end

function grid_points(dims::NTuple{3, Int}, origin::NTuple{3, <:Real}, spacing::NTuple{3, <:Real})
    points = Vector{SVector{3, Float64}}(undef, prod(dims))
    idx = 1
    for k in 1:dims[3], j in 1:dims[2], i in 1:dims[1]
        points[idx] = SVector{3, Float64}(
            Float64(origin[1]) + (i - 0.5) * Float64(spacing[1]),
            Float64(origin[2]) + (j - 0.5) * Float64(spacing[2]),
            Float64(origin[3]) + (k - 0.5) * Float64(spacing[3]),
        )
        idx += 1
    end
    return points
end

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