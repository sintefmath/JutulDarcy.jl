export WellDomain, MultiSegmentWell
export TotalMassFlux, PotentialDropBalanceWell, SegmentWellBoreFrictionHB

export InjectorControl, ProducerControl, SinglePhaseRateTarget, BottomHolePressureTarget

export Perforations
export MixedWellSegmentFlow
export segment_pressure_drop


abstract type WellPotentialFlowDiscretization <: PotentialFlowDiscretization end

"""
Two point approximation with flux for wells
"""
struct MixedWellSegmentFlow <: WellPotentialFlowDiscretization end

struct SurfaceWellConditions{T, R} <: ScalarVariable
    storage::T
    separator_conditions::Vector{NamedTuple{(:p, :T), Tuple{R, R}}}
    separator_targets::Vector{Tuple{Int, Int}}
    function SurfaceWellConditions(S::T, c, t, R::DataType = Float64) where T<:JutulStorage
        new{T, R}(S, c, t)
    end
end

function SurfaceWellConditions(sys::JutulSystem; kwarg...)
    s = JutulStorage()
    cond = [default_surface_cond()]
    targets = [(0, 0)]
    return SurfaceWellConditions(s, cond, targets)
end

include("separator.jl")

# Total velocity in each well segment
struct TotalMassFlux{R} <: ScalarVariable
    scale::Union{Nothing, R}
    max_abs::Union{Nothing, R}
    max_rel::Union{Nothing, R}
    function TotalMassFlux(;scale = 3600*24, max_abs = nothing, max_rel = nothing)
        new{Float64}(scale, max_abs, max_rel)
    end
end

relative_increment_limit(tmf::TotalMassFlux) = tmf.max_rel
absolute_increment_limit(tmf::TotalMassFlux) = tmf.max_abs

associated_entity(::TotalMassFlux) = Faces()
Jutul.variable_scale(t::TotalMassFlux) = t.scale

default_surface_cond() = (p = 101325.0, T = 288.15) # Pa and deg. K from ISO 13443:1996 for natural gas

function common_well_setup(nr; dz = nothing, WI = nothing, gravity = gravity_constant)
    if isnothing(dz)
        @warn "dz not provided for well. Assuming no gravity."
        gdz = zeros(nr)
    else
        @assert length(dz) == nr  "Must have one connection drop dz per perforated cell"
        gdz = dz*gravity
    end
    if isnothing(WI)
        @warn "No well indices provided. Using 1e-12."
        WI = repeat(1e-12, nr)
    else
        @assert length(WI) == nr  "Must have one well index per perforated cell ($(length(WI)) well indices, $nr reservoir cells))"
    end
    return (WI, gdz)
end

export setup_well, setup_vertical_well
"""
    setup_well(D::DataDomain, reservoir_cells; skin = 0.0, Kh = nothing, radius = 0.1, dir = :z)

Set up a well in `reservoir_cells` with given skin factor and radius. The order
of cells matter as it is treated as a trajectory.
"""
function setup_well(D::DataDomain, reservoir_cells; cell_centers = D[:cell_centroids], kwarg...)
    K = D[:permeability]
    g = physical_representation(D)
    return setup_well(g, K, reservoir_cells; cell_centers = cell_centers, kwarg...)
end

function setup_well(g, K, reservoir_cells::AbstractVector;
                                        reference_depth = nothing, 
                                        cell_centers = nothing,
                                        skin = 0.0,
                                        Kh = nothing,
                                        radius = 0.1,
                                        simple_well = false,
                                        WI = missing,
                                        dir = :z,
                                        kwarg...)
    n = length(reservoir_cells)
    # Make sure these are cell indices
    reservoir_cells = map(i -> cell_index(g, i), reservoir_cells)
    if isnothing(cell_centers)
        geometry = tpfv_geometry(g)
        cell_centers = geometry.cell_centroids
    end
    centers = cell_centers[:, reservoir_cells]

    if isnothing(reference_depth)
        reference_depth = centers[3, 1]
    end
    volumes = zeros(n)
    WI_computed = zeros(n)
    dz = zeros(n)

    function get_entry(x::AbstractVector, i)
        return x[i]
    end
    function get_entry(x, i)
        return x
    end
    for (i, c) in enumerate(reservoir_cells)
        if K isa AbstractVector
            k_i = K[c]
        else
            k_i = K[:, c]
        end
        WI_i = get_entry(WI, i)
        Kh_i = get_entry(Kh, i)
        r_i = get_entry(radius, i)
        s_i = get_entry(skin, i)
        if ismissing(WI_i) || isnan(WI_i)
            WI_i = compute_peaceman_index(g, k_i, r_i, c, skin = s_i, Kh = Kh_i)
        end
        WI_i::AbstractFloat
        WI_computed[i] = WI_i
        center = vec(centers[:, i])
        dz[i] = center[3] - reference_depth
        if dir isa Symbol
            d = dir
        else
            d = dir[i]
        end
        Δ = cell_dims(g, c)
        d_index = findfirst(isequal(d), [:x, :y, :z])
        h = Δ[d_index]
        volumes[i] = h*π*r_i^2
    end
    if simple_well
        W = SimpleWell(reservoir_cells, WI = WI_computed, dz = dz, reference_depth = reference_depth, kwarg...)
    else
        # Depth differences are taken care of via centers.
        dz *= 0.0
        W = MultiSegmentWell(reservoir_cells, volumes, centers; WI = WI_computed, dz = dz, reference_depth = reference_depth, kwarg...)
    end
    return W
end

function Jutul.plot_primitives(mesh::MultiSegmentWell, plot_type; kwarg...)
    # By default, no plotting is supported
    if plot_type == :lines
        centers = mesh.centers
        if ndims(centers) == 3
            # Some bug somewhere in MRST parser.
            centers = dropdims(centers, dims = 3)
        end
        pts = collect(centers')
        top = [pts[1, 1] pts[1, 2] pts[1, 3] - 10.0]
        pts = vcat(top, pts)
        @. pts[:, 3] *= -1

        function cell_mapper(x::AbstractVector)
            return vcat(x[1], x)
        end

        function cell_mapper(x::AbstractMatrix)
            return hcat(x[:, 1], x)
        end
        mapper = (Cells = x -> cell_mapper(x), )
        out = (points = pts, mapper = mapper, top_text = String(mesh.name), marker_size = 20)
    else
        out = nothing
    end
    return out
end

"""
    setup_vertical_well(D::DataDomain, i, j; <kwarg>)

Set up a vertical well with a `DataDomain` input that represents the porous
medium / reservoir where the wells it to be placed.
"""
function setup_vertical_well(D::DataDomain, i, j; cell_centers = D[:cell_centroids], kwarg...)
    K = D[:permeability]
    g = physical_representation(D)
    return setup_vertical_well(g, K, i, j; cell_centers = cell_centers, kwarg...)
end

"""
    setup_vertical_well(g, K, i, j; heel = 1, toe = grid_dims_ijk(g)[3], kwarg...)

Set up a vertical well for given grid `g` and permeability `K` at logical
indices `i, j` perforating all cells starting at k-logical index `heel` to
`toe`.
"""
function setup_vertical_well(g, K, i, j; heel = 1, toe = grid_dims_ijk(g)[3], kwarg...)
    @assert heel <= toe
    @assert heel > 0
    @assert toe > 0
    k_range = heel:toe
    n = length(k_range)
    @assert n > 0
    reservoir_cells = zeros(Int64, n)
    for (ix, k) in enumerate(k_range)
        reservoir_cells[ix] = cell_index(g, (i, j, k))
    end
    return setup_well(g, K, reservoir_cells; kwarg...)
end

include("mswells_equations.jl")

function update_before_step_well!(well_state, well_model, res_state, res_model, ctrl)

end

function domain_fluid_volume(grid::WellDomain)
    return grid.volumes
end

function domain_fluid_volume(grid::SimpleWell)
    return [grid.volume]
end

# Well segments

function get_neighborship(::SimpleWell)
    # No interior connections.
    return zeros(Int64, 2, 0)
end

function number_of_cells(W::SimpleWell)
    return 1
end

function number_of_cells(W::WellDomain)
    return length(W.volumes)
end

function declare_entities(W::WellDomain)
    np = length(W.perforations.self)
    c = (entity = Cells(),         count = number_of_cells(W))
    f = (entity = Faces(),         count = number_of_faces(W))
    p = (entity = Perforations(),  count = np)
    return [c, f, p]
end

function Jutul.select_secondary_variables!(S, D::WellDomain, model)
    S[:SurfaceWellConditions] = SurfaceWellConditions(model.system)
end

Base.@propagate_inbounds function multisegment_well_perforation_flux!(out, sys::Union{ImmiscibleSystem, SinglePhaseSystem}, state_res, state_well, rhoS, conn)
    rc = conn.reservoir
    wc = conn.well
    # Reservoir quantities
    ρ = state_res.PhaseMassDensities
    # Extra mobility needed
    kr = state_res.RelativePermeabilities
    μ = state_res.PhaseViscosities
    nph = size(ρ, 1)
    # Well quantities
    ρ_w = state_well.PhaseMassDensities
    # Saturation instead of mobility - use total mobility form
    s_w = state_well.Saturations
    λ_t = 0
    for ph in 1:nph
        λ_t += kr[ph, rc]/μ[ph, rc]
    end
    for ph in 1:nph
        # ψ is pressure difference from reservoir to well. If it is negative, we are injecting into the reservoir.
        ψ = perforation_phase_potential_difference(conn, state_res, state_well, ph)
        if ψ < 0
            # Injection
            out[ph] = s_w[ph, wc]*ρ_w[ph, wc]*ψ*λ_t
        else
            # Production
            λ = kr[ph, rc]/μ[ph, rc]
            out[ph] = λ*ρ[ph, rc]*ψ
        end
    end
    return out
end

Base.@propagate_inbounds function simple_well_perforation_flux!(out, sys::Union{ImmiscibleSystem, SinglePhaseSystem}, state_res, state_well, rhoS, conn)
    rc = conn.reservoir
    # Reservoir quantities
    ρ = state_res.PhaseMassDensities
    # Extra mobility needed
    kr = state_res.RelativePermeabilities
    μ = state_res.PhaseViscosities
    nph = size(ρ, 1)
    ρλ_t = 0
    for ph in 1:nph
        ρλ_t += ρ[ph, rc]*kr[ph, rc]/μ[ph, rc]
    end
    X = state_well.MassFractions
    for ph in 1:nph
        # ψ is pressure difference from reservoir to well. If it is negative, we are injecting into the reservoir.
        ψ = perforation_phase_potential_difference(conn, state_res, state_well, ph)
        if ψ < 0
            # Injection
            out[ph] = X[ph]*ψ*ρλ_t
        else
            # Production
            λ = kr[ph, rc]/μ[ph, rc]
            out[ph] = λ*ρ[ph, rc]*ψ
        end
    end
    return out
end

include("mswells.jl")
include("stdwells.jl")

include("stdwells_equations.jl")

# Some utilities
function mix_by_mass(masses, total::T, values) where T
    v = zero(T)
    @inbounds for i in eachindex(masses)
        v += masses[i]*values[i]
    end
    return v/total
end

function mix_by_saturations(s, values)
    v = zero(eltype(s))
    @inbounds for i in eachindex(s)
        v += s[i]*values[i]
    end
    return v
end

function mix_by_saturations(s::Real, values)
    return s*values[]
end

function setup_forces(model::SimulationModel{D, S}; mask = nothing) where {D <: DiscretizedDomain{G}, S<:MultiPhaseSystem} where G<:WellDomain
    mask::Union{Nothing, PerforationMask}
    return (mask = mask,)
end

function apply_perforation_mask!(storage::NamedTuple, mask::AbstractVector)
    function mask_row!(M::AbstractMatrix, m, ix)
        for i in axes(M, 1)
            M[i, ix] *= m
        end
    end
    function mask_row!(M::AbstractVector, m, ix)
        M[ix] *= m
    end
    for (k, s) in pairs(storage)
        if k == :numeric
            continue
        end
        v = s.entries
        for i in 1:Jutul.number_of_entities(s)
            mask_value = mask[i]
            for j in Jutul.vrange(s, i)
                mask_row!(v, mask_value, j)
            end
        end
    end
end

function flash_wellstream_at_surface(well_model, well_state, rhoS)
    fsys = flow_system(well_model.system)
    return flash_wellstream_at_surface(well_model, fsys, well_state, rhoS)
end

function flash_wellstream_at_surface(well_model, system::ImmiscibleSystem, well_state, rhoS)
    vol = well_state.TotalMasses[:, 1]./rhoS
    volfrac = vol./sum(vol)
    return (rhoS, volfrac)
end

function flash_wellstream_at_surface(well_model, system::SinglePhaseSystem, well_state, rhoS)
    return (rhoS, [1.0])
end

function surface_density_and_volume_fractions(state)
    x = only(state.SurfaceWellConditions)
    return (x.density, x.volume_fractions)
end
