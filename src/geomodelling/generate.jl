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
    box_origin::NTuple{3, <:Real} = (0.0, 0.0, 0.0),
    porosity_mean::Real = 0.20,
    porosity_std::Real = 0.05,
    porosity_bounds::Tuple{<:Real, <:Real} = (0.05, 0.95),
    permeability_bounds::Tuple{<:Real, <:Real} = (1e-20, Inf),
    perm_from_poro = kozeny_carman_permeability,
    porosity_process = _standard_normal_generator,
)
    nrealizations > 0 || throw(ArgumentError("nrealizations must be positive."))
    all(>(0), dims) || throw(ArgumentError("All grid dimensions must be positive."))
    all(>(0), box_lengths) || throw(ArgumentError("All box lengths must be positive."))

    phimin, phimax = porosity_bounds
    0.0 <= phimin < phimax < 1.0 || throw(ArgumentError("porosity_bounds must satisfy 0 <= min < max < 1."))

    if !isnothing(seed)
        Random.seed!(seed)
    end

    spacing = ntuple(i -> Float64(box_lengths[i]) / dims[i], 3)
    volumes = fill(Float64(prod(spacing)), dims...)
    points = grid_points(dims, box_origin, spacing)

    realizations = NamedTuple[]
    for _ in 1:nrealizations
        zporo = porosity_process(dims...)
        porosity = clamp.(porosity_mean .+ porosity_std .* zporo, phimin, phimax)
        permeability = perm_from_poro(porosity; bounds = permeability_bounds)
        push!(realizations, (
            porosity = porosity,
            permeability = permeability,
            points = points,
            volumes = volumes,
        ))
    end
    return nrealizations == 1 ? only(realizations) : realizations
end

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