export WellGrid, MultiSegmentWell
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

abstract type WellGrid <: PorousMediumGrid
    # Wells are not porous themselves per se, but they are discretizing
    # part of a porous medium.
end

function Base.show(io::IO, w::WellGrid)
    if w isa SimpleWell
        nseg = 0
    else
        nseg = size(w.neighborship, 2)
    end
    print(io, "$(typeof(w)) [$(w.name)] ($(length(w.volumes)) nodes, $(nseg) segments, $(length(w.perforations.reservoir)) perforations)")
end

# Total velocity in each well segment
struct TotalMassFlux <: ScalarVariable
    scale
    max_abs
    max_rel
    function TotalMassFlux(;scale = 3600*24, max_abs = nothing, max_rel = nothing)
        new(scale, max_abs, max_rel)
    end
end

relative_increment_limit(tmf::TotalMassFlux) = tmf.max_rel
absolute_increment_limit(tmf::TotalMassFlux) = tmf.max_abs

associated_entity(::TotalMassFlux) = Faces()
Jutul.variable_scale(t::TotalMassFlux) = t.scale

default_surface_cond() = (p = 101325.0, T = 288.15) # Pa and deg. K from ISO 13443:1996 for natural gas

struct SimpleWell <: WellGrid
    volumes
    perforations
    surface
    reservoir_symbol
    name::Symbol     # Symbol that names the well
    function SimpleWell(reservoir_cells; name = :Well, reference_depth = 0, volume = 1e-3, reservoir_symbol = :Reservoir, surface_conditions = default_surface_cond(), kwarg...)
        nr = length(reservoir_cells)

        WI, gdz = common_well_setup(nr; kwarg...)
        perf = (self = ones(Int64, nr), reservoir = vec(reservoir_cells), WI = WI, gdz = gdz)
        new([volume], perf, surface_conditions, reservoir_symbol, name)
    end
end

struct MultiSegmentWell <: WellGrid
    volumes          # One per cell
    perforations     # (self -> local cells, reservoir -> reservoir cells, WI -> connection factor)
    neighborship     # Well cell connectivity
    top              # "Top" node where scalar well quantities live
    centers          # Coordinate centers of nodes
    surface          # p, T at surface
    reservoir_symbol # Symbol of the reservoir the well is coupled to
    name::Symbol     # Symbol that names the well
    segment_models   # Segment pressure drop model
    function MultiSegmentWell(reservoir_cells, volumes::AbstractVector, centers;
                                                        N = nothing,
                                                        name = :Well,
                                                        perforation_cells = nothing,
                                                        segment_models = nothing,
                                                        reference_depth = 0,
                                                        dz = nothing,
                                                        surface_conditions = default_surface_cond(),
                                                        accumulator_volume = mean(volumes),
                                                        reservoir_symbol = :Reservoir, kwarg...)
        nv = length(volumes)
        nc = nv + 1
        reservoir_cells = vec(reservoir_cells)
        nr = length(reservoir_cells)
        if isnothing(N)
            @debug "No connectivity. Assuming nicely ordered linear well."
            N = vcat((1:nv)', (2:nc)')
        elseif maximum(N) == nv
            N = vcat([1, 2], N+1)
        end
        nseg = size(N, 2)
        @assert size(N, 1) == 2
        @assert size(centers, 1) == 3
        volumes = vcat([accumulator_volume], volumes)
        ext_centers = hcat([centers[1:2, 1]..., reference_depth], centers)
        @assert length(volumes) == size(ext_centers, 2)
        if !isnothing(reservoir_cells) && isnothing(perforation_cells)
            @assert length(reservoir_cells) == nv "If no perforation cells are given, we must 1->1 correspondence between well volumes and reservoir cells."
            perforation_cells = collect(2:nc)
        end
        perforation_cells = vec(perforation_cells)

        if isnothing(segment_models)
            Δp = SegmentWellBoreFrictionHB(1.0, 1e-4, 0.1)
            segment_models = repeat([Δp], nseg)
        else
            segment_models::AbstractVector
            @assert length(segment_models) == nseg
        end
        if isnothing(dz)
            dz = centers[3, :] - reference_depth
        end
        @assert length(perforation_cells) == nr
        WI, gdz = common_well_setup(nr; dz = dz, kwarg...)
        perf = (self = perforation_cells, reservoir = reservoir_cells, WI = WI, gdz = gdz)
        accumulator = (reference_depth = reference_depth, )
        new(volumes, perf, N, accumulator, ext_centers, surface_conditions, reservoir_symbol, name, segment_models)
    end
end

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
        @assert length(WI) == nr  "Must have one well index per perforated cell"
    end
    return (WI, gdz)
end

export setup_well, setup_vertical_well
function setup_well(g, K, reservoir_cells::AbstractVector;
                                        reference_depth = nothing, 
                                        geometry = tpfv_geometry(g),
                                        skin = 0.0,
                                        Kh = nothing,
                                        radius = 0.1,
                                        simple_well = false,
                                        dir = :z,
                                        kwarg...)
    n = length(reservoir_cells)
    # Make sure these are cell indices
    reservoir_cells = map(i -> cell_index(g, i), reservoir_cells)
    centers = geometry.cell_centroids[:, reservoir_cells]
    if isnothing(reference_depth)
        reference_depth = centers[3, 1]
    end
    volumes = zeros(n)
    WI = zeros(n)
    dz = zeros(n)
    for (i, c) in enumerate(reservoir_cells)
        if K isa AbstractVector
            k = K[c]
        else
            k = K[:, c]
        end
        WI[i] = compute_peaceman_index(g, k, radius, c, skin = skin, Kh = Kh)
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
        volumes[i] = h*π*radius^2
    end
    if simple_well
        W = SimpleWell(reservoir_cells, WI = WI, dz = dz, reference_depth = reference_depth, kwarg...)
    else
        W = MultiSegmentWell(reservoir_cells, volumes, centers; WI = WI, dz = dz, reference_depth = reference_depth, kwarg...)
    end
    return W
end

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

"""
Hagedorn and Brown well bore friction model for a segment.
"""
struct SegmentWellBoreFrictionHB{R}
    L::R
    roughness::R
    D_outer::R
    D_inner::R
    assume_turbulent::Bool
    laminar_limit::R
    turbulent_limit::R
    function SegmentWellBoreFrictionHB(L, roughness, D_outer; D_inner = 0, assume_turbulent = false, laminar_limit = 2000.0, turbulent_limit = 2000.0)
        new{typeof(L)}(L, roughness, D_outer, D_inner, assume_turbulent, laminar_limit, turbulent_limit)
    end
end

function is_turbulent_flow(f::SegmentWellBoreFrictionHB, Re)
    return f.assume_turbulent || Re >= f.turbulent_limit
end

function is_laminar_flow(f::SegmentWellBoreFrictionHB, Re)
    return !f.assume_turbulent && Re <= f.laminar_limit
end

function segment_pressure_drop(f::SegmentWellBoreFrictionHB, v, ρ, μ)
    D⁰, Dⁱ = f.D_outer, f.D_inner
    R, L = f.roughness, f.L
    ΔD = D⁰-Dⁱ
    s = v > 0.0 ? 1.0 : -1.0
    e = eps(Float64)
    v = s*max(abs(v), e)
    # Scaling - assuming input is total mass rate
    v = v/(π*ρ*((D⁰/2)^2 - (Dⁱ/2)^2))
    Re = abs(v*ρ*ΔD)/μ
    # Friction model - empirical relationship
    Re_l, Re_t = f.laminar_limit, f.turbulent_limit
    if is_laminar_flow(f, Re)
        f = 16/Re_l
    else
        # Either turbulent or intermediate flow regime. We need turbulent value either way.
        f_t = (-3.6*log(6.9/Re +(R/(3.7*D⁰))^(10/9))/log(10))^(-2)
        if is_turbulent_flow(f, Re)
            # Turbulent flow
            f = f_t
        else
            # Intermediate regime - interpolation
            f_l = 16/Re_l
            Δf = f_t - f_l
            ΔRe = Re_t - Re_l
            f = f_l + (Δf / ΔRe)*(Re - Re_l)
        end
    end
    Δp = -(2*s*L/ΔD)*(f*ρ*v^2)
    return Δp
end


struct PotentialDropBalanceWell{T} <: JutulEquation
    flow_discretization::T
end

associated_entity(::PotentialDropBalanceWell) = Faces()

include("well_equations.jl")

function convergence_criterion(model, storage, eq::PotentialDropBalanceWell, eq_s, r; dt = 1)
    e = [norm(r, Inf)/1e5] # Given as pressure - scale by 1 bar
    R = Dict("AbsMax" => (errors = e, names = "R"))
    return R
end


function fluid_volume(grid::WellGrid)
    return grid.volumes
end

# Well segments
"""
Perforations are connections from well cells to reservoir vcells
"""
struct Perforations <: JutulUnit end

function get_neighborship(::SimpleWell)
    # No interior connections.
    return zeros(Int64, 2, 0)
end

function get_neighborship(W::MultiSegmentWell)
    return W.neighborship
end

function number_of_cells(W::WellGrid)
    length(W.volumes)
end

function declare_entities(W::WellGrid)
    c = (entity = Cells(),         count = number_of_cells(W))
    f = (entity = Faces(),         count = number_of_faces(W))
    p = (entity = Perforations(),  count = length(W.perforations.self))
    return [c, f, p]
end

Base.@propagate_inbounds function well_perforation_flux!(out, sys::Union{ImmiscibleSystem, SinglePhaseSystem}, state_res, state_well, rhoS, dp, rc, wc)
    # Reservoir quantities
    ρ = state_res.PhaseMassDensities
    # Extra mobility needed
    kr = state_res.RelativePermeabilities
    μ = state_res.PhaseViscosities
    # Well quantities
    ρ_w = state_well.PhaseMassDensities
    # Saturation instead of mobility - use total mobility form
    s_w = state_well.Saturations
    nph = size(s_w, 1)
    # dp is pressure difference from reservoir to well. If it is negative, we are injecting into the reservoir.
    if dp < 0
        # Injection
        λ_t = 0
        for ph in 1:nph
            λ_t += kr[ph, rc]/μ[ph, rc]
        end
        Q = λ_t*dp
        for ph in 1:nph
            out[ph] = s_w[ph, wc]*ρ_w[ph, wc]*Q
        end
    else
        # Production
        for ph in 1:nph
            λ = kr[ph, rc]/μ[ph, rc]
            out[ph] = λ*ρ[ph, rc]*dp
        end
    end
end

Base.@propagate_inbounds function well_perforation_flux!(out, sys::CompositionalSystem, state_res, state_well, rhoS, dp, rc, wc)
    μ = state_res.PhaseViscosities
    kr = state_res.RelativePermeabilities
    ρ = state_res.PhaseMassDensities
    X = state_res.LiquidMassFractions
    Y = state_res.VaporMassFractions

    ρ_w = state_well.PhaseMassDensities
    s_w = state_well.Saturations
    X_w = state_well.LiquidMassFractions
    Y_w = state_well.VaporMassFractions

    nc = size(X, 1)
    nph = size(μ, 1)
    has_water = nph == 3
    phase_ix = phase_indices(sys)

    if has_water
        A, L, V = phase_ix
    else
        L, V = phase_ix
    end
    # dp is pressure difference from reservoir to well. If it is negative, we are injecting into the reservoir.
    mob(ph) = kr[ph, rc]/μ[ph, rc]
    if dp < 0
        # Injection
        λ_t = 0
        for ph in 1:nph
            λ_t += mob(ph)
        end
        Q = λ_t*dp
        phase_mass_flux = ph -> ρ_w[ph, wc]*s_w[ph, wc]*Q
        Q_l = phase_mass_flux(L)
        Q_v = phase_mass_flux(V)
        X_upw = X_w
        Y_upw = Y_w
        cell_upw = wc
    else
        phase_mass_flux = ph -> ρ[ph, rc]*mob(ph)*dp
        X_upw = X
        Y_upw = Y
        cell_upw = rc
    end
    Q_l = phase_mass_flux(L)
    Q_v = phase_mass_flux(V)

    @inbounds for c in 1:nc
        out[c] = Q_l*X_upw[c, cell_upw] + Q_v*Y_upw[c, cell_upw]
    end
    if has_water
        out[nc+1] = phase_mass_flux(A)
    end
end

# Selection of primary variables
function select_primary_variables_domain!(S, domain::DiscretizedDomain{G}, system, formulation) where {G<:MultiSegmentWell}
    S[:TotalMassFlux] = TotalMassFlux()
end

function select_equations_domain!(eqs, domain::DiscretizedDomain{G}, system, arg...) where {G<:MultiSegmentWell}
    eqs[:potential_balance] = PotentialDropBalanceWell(domain.discretizations.mass_flow)
end

function minimum_output_variables(domain::DiscretizedDomain{G}, system::CompositionalSystem, formulation, primary_variables, secondary_variables) where {G<:WellGrid}
    vars = minimum_output_variables(system, primary_variables)
    push!(vars, :PhaseMassDensities)
    push!(vars, :Saturations)
    return vars
end

# Some utilities
function mix_by_mass(masses, total, values)
    v = 0
    @inbounds for i in eachindex(masses)
        v += masses[i]*values[i]
    end
    return v/total
end

function mix_by_saturations(s, values)
    v = 0
    @inbounds for i in eachindex(s)
        v += s[i]*values[i]
    end
    return v
end

function mix_by_saturations(s::Real, values)
    return s*values[]
end

function setup_forces(model::SimulationModel{D, S}; mask = nothing) where {D <: DiscretizedDomain{G}, S<:MultiPhaseSystem} where G<:WellGrid
    mask::Union{Nothing, PerforationMask}
    return (mask = mask,)
end

function apply_mask!(v::AbstractMatrix, pm::PerforationMask)
    for (i, mask) in enumerate(pm.values)
        for r in 1:size(v, 1)
            v[r, i] *= mask
        end
    end
end

function apply_mask!(v::AbstractVector, pm::PerforationMask)
    for (i, mask) in enumerate(pm.values)
        v[i] *= mask
    end
end

function apply_mask!(ct::InjectiveCrossTerm, pm::PerforationMask)
    apply_mask!(ct.crossterm_target, pm)
    apply_mask!(ct.crossterm_source, pm)
end

function Jutul.apply_force_to_cross_term_target!(ct::InjectiveCrossTerm, storage_t, storage_s, model_t, model_s, source, target, force::PerforationMask, dt, time)
    if source == :Reservoir || target == :Reservoir
        apply_mask!(ct, force)
    end
end

function Jutul.apply_force_to_cross_term_source!(ct::InjectiveCrossTerm, storage_t, storage_s, model_t, model_s, source, target, force::PerforationMask, dt, time)
    if source == :Reservoir || target == :Reservoir
        apply_mask!(ct, force)
    end
end
