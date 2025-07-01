

abstract type WellPotentialFlowDiscretization <: PotentialFlowDiscretization end

"""
Two point approximation with flux for wells
"""
struct MixedWellSegmentFlow <: WellPotentialFlowDiscretization end

include("separator.jl")
include("well_results.jl")

"""
    TotalMassFlux(scale = si_unit(:day), max_abs = nothing, max_rel = nothing)

Variable normally used as primary variable. Represents the total mass flux going
through a face. The typical usage is the mass flow through a segment of a
[`MultiSegmentWell`](@ref).

Note that the flow direction can often switch signs over a segment during a
complex simulation. Setting `max_rel` to something other than `nothing` can
therefore lead to severe convergence issues in the case of flow reversal.

# Fields (as keyword arguments)

$FIELDS
"""
struct TotalMassFlux{R} <: ScalarVariable
    "Scaling for variable"
    scale::Union{Nothing, R}
    "Max absolute change between Newton iterations"
    max_abs::Union{Nothing, R}
    "Maximum relative change between Newton iterations"
    max_rel::Union{Nothing, R}
    function TotalMassFlux(; scale = si_unit(:day), max_abs = nothing, max_rel = nothing)
        new{Float64}(scale, max_abs, max_rel)
    end
end

relative_increment_limit(tmf::TotalMassFlux) = tmf.max_rel
absolute_increment_limit(tmf::TotalMassFlux) = tmf.max_abs

associated_entity(::TotalMassFlux) = Faces()
Jutul.variable_scale(t::TotalMassFlux) = t.scale

default_surface_cond() = (p = 101325.0, T = 288.15) # Pa and deg. K from ISO 13443:1996 for natural gas

function common_well_setup(nr; dz = nothing, WI = nothing, WIth = nothing, gravity = gravity_constant)
    if isnothing(dz)
        @warn "dz not provided for well. Assuming no gravity."
        gdz = zeros(nr)
    else
        @assert length(dz) == nr  "Must have one connection drop dz per perforated cell"
        gdz = dz*gravity
    end
    if isnothing(WI)
        @warn "No well indices provided. Using 1e-12."
        WI = fill(1e-12, nr)
    else
        @assert length(WI) == nr  "Must have one well index per perforated cell ($(length(WI)) well indices, $nr reservoir cells))"
    end
    if isnothing(WIth)
        @warn "No thermal well indices provided. Using 1."
        WIth = fill(1.0, nr)
    else
        @assert length(WIth) == nr  "Must have one thermal well index per perforated cell ($(length(WIth)) thermal well indices, $nr reservoir cells))"
    end
    return (WI, WIth, gdz)
end

"""
    setup_well(D::DataDomain, reservoir_cells; skin = 0.0, Kh = nothing, radius = 0.1, dir = :z, name = :Well)
    w = setup_well(D, 1, name = :MyWell)         # Cell 1 in the grid
    w = setup_well(D, (2, 5, 1), name = :MyWell) # Position (2, 5, 1) in logically structured mesh
    w2 = setup_well(D, [1, 2, 3], name = :MyOtherWell)


Set up a well in `reservoir_cells` with given skin factor and radius. The order
of cells matter as it is treated as a trajectory.

The `name` keyword argument can be left defaulted if your model will only have a
single well (named `:Well`). It is highly recommended to provide this whenever a
well is set up.

`reservoir_cells` can be one of the following: A Vector of cells, a single cell,
a Vector of `(I, J, K)` Tuples or a single Tuple of the same type.
"""
function setup_well(D::DataDomain, reservoir_cells; cell_centers = D[:cell_centroids], kwarg...)
    # Get permeability
    K = D[:permeability]
    # Compute effective thermal conductivity
    Λ_f = D[:fluid_thermal_conductivity]
    Λ_r = D[:rock_thermal_conductivity]
    ϕ = D[:porosity]
    Λ_r = vec(Λ_r)
    ϕ = vec(ϕ)
    if size(Λ_f, 1) == 1 || Λ_f isa Vector
        Λ_f = vec(Λ_f)
        Λ = ϕ.*Λ_f + (1.0 .- ϕ).*Λ_r
    else
        # TODO: This is a bit of a hack. We should really have a proper way to
        # do this inside the equations for multiphase flow.
        Λ = Λ_r
    end
    # Get grid
    g = physical_representation(D)
    return setup_well(g, K, reservoir_cells; thermal_conductivity = Λ, cell_centers = cell_centers, kwarg...)
end

function setup_well(g, K, reservoir_cells::AbstractVector;
        reference_depth = nothing,
        cell_centers = nothing,
        skin = 0.0,
        Kh = nothing,
        radius = 0.1,
        accumulator_volume = missing,
        simple_well = true,
        simple_well_regularization = 1.0,
        WI = missing,
        WIth = missing,
        thermal_conductivity = missing,
        thermal_index_args = NamedTuple(),
        dir = :z,
        kwarg...
    )
    T = promote_type(eltype(K), eltype(skin), eltype(radius), typeof(simple_well_regularization))
    T = promote_type(T, Jutul.float_type(g))
    n = length(reservoir_cells)
    # Make sure these are cell indices
    reservoir_cells = map(i -> cell_index(g, i), reservoir_cells)
    if isnothing(cell_centers)
        geometry = tpfv_geometry(g)
        cell_centers = geometry.cell_centroids
    end
    centers = cell_centers[:, reservoir_cells]

    if isnothing(reference_depth)
        if size(centers, 1) == 2
            reference_depth = 0.0
        else
            reference_depth = centers[3, 1]
        end
    end
    volumes = zeros(T, n)
    WI_computed = zeros(T, n)
    WIth_computed = zeros(T, n)
    Λ = thermal_conductivity
    dz = zeros(T, n)

    function get_entry(x::AbstractVector, i)
        return x[i]
    end
    function get_entry(x, i)
        return x
    end
    segment_radius = Float64[]
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
            WI_i = compute_peaceman_index(g, k_i, r_i, c, skin = s_i, Kh = Kh_i, dir = dir)
        end
        WI_computed[i] = WI_i
        WIth_i = 0.0
        if !ismissing(Λ)
            WIth_i = get_entry(WIth, i)
            if Λ isa AbstractVector
                Λ_i = Λ[c]
            else
                Λ_i = Λ[:, c]
            end
            if ismissing(WIth_i) || isnan(WIth_i)
                WIth_i = compute_well_thermal_index(g, Λ_i, r_i, c; 
                    dir = dir, thermal_index_args...)
            end
        end
        push!(segment_radius, r_i)
        WIth_computed[i] = WIth_i
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
        if ismissing(accumulator_volume)
            accumulator_volume = simple_well_regularization*maximum(volumes)
        end
        W = SimpleWell(reservoir_cells; WI = WI_computed, WIth = WIth_computed, volume = accumulator_volume, dz = dz, reference_depth = reference_depth, kwarg...)
    else
        # Depth differences are taken care of via centers.
        dz *= 0.0
        W = MultiSegmentWell(reservoir_cells, volumes, centers; 
            WI = WI_computed, WIth = WIth_computed, dz = dz, 
            reference_depth = reference_depth, segment_radius = segment_radius, 
            kwarg...)
    end
    return W
end

function setup_well(g, K, reservoir_cell::Union{Int, Tuple, NamedTuple}; kwarg...)
    return setup_well(g, K, [reservoir_cell]; kwarg...)
end

function setup_well_from_trajectory(D::DataDomain, traj; traj_arg = NamedTuple(), kwarg...)
    G = D |> physical_representation |> UnstructuredMesh
    cells = Jutul.find_enclosing_cells(G, traj; traj_arg...)
    return setup_well(D, cells; kwarg...)
end

function map_well_nodes_to_reservoir_cells(w::MultiSegmentWell, reservoir::Union{DataDomain, Missing} = missing)
    # TODO: Try to more or less match it up cell by cell. Could be
    # improved...
    c = zeros(Int, length(w.volumes))
    c[w.perforations.self] .= w.perforations.reservoir
    for i in 2:length(c)
        if c[i] == 0
            c[i] = c[i-1]
        end
    end
    for i in (length(c)-1):-1:1
        if c[i] == 0
            c[i] = c[i+1]
        end
    end
    @assert all(x -> x > 0, c)
    return c
end

function map_well_nodes_to_reservoir_cells(w::SimpleWell, reservoir::Union{DataDomain, Missing} = missing)
    return [w.perforations.reservoir[1]]
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
    setup_vertical_well(D::DataDomain, i, j; name = :MyWell, <kwarg>)

Set up a vertical well with a [`DataDomain`](@ref) input that represents the porous
medium / reservoir where the wells it to be placed. See [`SimpleWell`](@ref),
[`MultiSegmentWell`](@ref) and [`setup_well`](@ref) for more details about possible keyword
arguments.
"""
function setup_vertical_well(D::DataDomain, i, j; cell_centers = D[:cell_centroids], kwarg...)
    # Get permeability
    K = D[:permeability]
    # Compute effective thermal conductivity
    Λ_f = D[:fluid_thermal_conductivity]
    Λ_r = D[:rock_thermal_conductivity]
    ϕ = D[:porosity]
    Λ = ϕ.*Λ_f + (1.0 .- ϕ).*Λ_r
    # Get grid
    g = physical_representation(D)
    return setup_vertical_well(g, K, i, j; thermal_conductivity = Λ, cell_centers = cell_centers, kwarg...)
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
    reservoir_cells = Int[]
    for (ix, k) in enumerate(k_range)
        cell_ix = cell_index(g, (i, j, k), throw = false)
        if isnothing(cell_ix)
            jutul_message("Well", "Cell ($i, $j, $k) not found in active set, skipping.", color = :yellow)
        else
            push!(reservoir_cells, cell_ix)
        end
    end
    length(reservoir_cells) > 0 || error("No cells found for well.")
    return setup_well(g, K, reservoir_cells; kwarg...)
end

include("mswells_equations.jl")

function update_before_step_well!(well_state, well_model, res_state, res_model, ctr, mask; kwarg...)

end

function domain_fluid_volume(grid::WellDomain)
    return grid.volumes.*grid.void_fraction
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
    if model.system isa MultiPhaseSystem || model.system isa CompositeSystem
        sys = flow_system(model.system)
        S[:SurfaceWellConditions] = SurfaceWellConditions(sys)
    end
end

Base.@propagate_inbounds function multisegment_well_perforation_flux!(out, sys::Union{ImmiscibleSystem, SinglePhaseSystem}, state_res, state_well, rhoS, conn)
    rc = conn.reservoir
    wc = conn.well
    nph = number_of_phases(sys)
    λ_t = sum(perforation_reservoir_mobilities(state_res, state_well, sys, rc, wc))
    for ph in 1:nph
        out[ph] = perforation_phase_mass_flux(λ_t, conn, state_res, state_well, ph)
    end
    return out
end

function perforation_phase_mass_flux(λ_t, conn, state_res, state_well, ph)
    # ψ is pressure difference from reservoir to well. If it is negative, we are injecting into the reservoir.
    ψ = perforation_phase_potential_difference(conn, state_res, state_well, ph)
    if ψ < 0
        wc = conn.well
        # Injection
        ρ_w = state_well.PhaseMassDensities
        s_w = state_well.Saturations
        q_ph = s_w[ph, wc]*ρ_w[ph, wc]*ψ*λ_t
    else
        rc = conn.reservoir
        # Production
        ρ = state_res.PhaseMassDensities
        λ = state_res.PhaseMobilities[ph, rc]
        q_ph = λ*ρ[ph, rc]*ψ
    end
    return q_ph
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

function apply_perforation_mask!(M::AbstractVector, mask::AbstractVector)
    for i in eachindex(mask)
        M[i] *= mask[i]
    end
    return M
end

function apply_perforation_mask!(M::AbstractMatrix, mask::AbstractVector)
    for j in eachindex(mask)
        for i in axes(M, 1)
            M[i, j] *= mask[j]
        end
    end
    return M
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

function flash_wellstream_at_surface(var, well_model, well_state, rhoS, cond = default_surface_cond())
    fsys = flow_system(well_model.system)
    return flash_wellstream_at_surface(var, well_model, fsys, well_state, rhoS)
end

function flash_wellstream_at_surface(var, well_model, system::ImmiscibleSystem, well_state, rhoS, cond = default_surface_cond())
    vol = well_state.TotalMasses[:, 1]./rhoS
    volfrac = vol./sum(vol)
    return (rhoS, volfrac)
end

function flash_wellstream_at_surface(var, well_model, system::SinglePhaseSystem, well_state, rhoS, cond = default_surface_cond())
    return (rhoS, [1.0])
end

function surface_density_and_volume_fractions(state)
    x = only(state.SurfaceWellConditions)
    return (x.density, x.volume_fractions)
end
