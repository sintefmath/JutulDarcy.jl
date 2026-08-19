function _standard_normal_generator(dims...)
    return randn(dims...)
end

"""
    generate_perm_poro(dims, box_lengths = (1.0, 1.0, 1.0); kwargs...)

Generate porosity and permeability fields from a latent Gaussian draw.

The default generator samples standard normal draws and applies the affine
scaling and translation implied by `porosity_mean` and `porosity_std` before
clipping to `porosity_bounds`. Permeability is then derived with the function
passed as `perm_from_poro`, which defaults to `kozeny_carman_permeability`.
"""
function generate_perm_poro(
    dims::NTuple{3, Int},
    box_lengths::NTuple{3, <:Real} = (1.0, 1.0, 1.0);
    nrealizations::Int = 1,
    seed = nothing,
    porosity_mean::Real = 0.20,
    porosity_std::Real = 0.05,
    porosity_bounds::Tuple{<:Real, <:Real} = (0.05, 0.95),
    permeability_bounds::Tuple{<:Real, <:Real} = (1e-20, Inf),
    perm_from_poro = kozeny_carman_permeability,
    generator = _standard_normal_generator,
)
    nrealizations > 0 || throw(ArgumentError("nrealizations must be positive."))
    all(>(0), dims) || throw(ArgumentError("All grid dimensions must be positive."))
    all(>(0), box_lengths) || throw(ArgumentError("All box lengths must be positive."))

    phimin, phimax = porosity_bounds
    0.0 <= phimin < phimax < 1.0 || throw(ArgumentError("porosity_bounds must satisfy 0 <= min < max < 1."))

    if !isnothing(seed)
        Random.seed!(seed)
    end

    realizations = NamedTuple[]
    for _ in 1:nrealizations
        zporo = generator(dims...)
        porosity = clamp.(porosity_mean .+ porosity_std .* zporo, phimin, phimax)
        permeability = perm_from_poro(porosity; bounds = permeability_bounds)
        push!(realizations, (porosity = porosity, permeability = permeability))
    end
    return nrealizations == 1 ? only(realizations) : realizations
end

function generate_perm_poro_realizations(dims::NTuple{3, Int}, box_lengths::NTuple{3, <:Real} = (1.0, 1.0, 1.0); kwargs...)
    return generate_perm_poro(dims, box_lengths; kwargs...)
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