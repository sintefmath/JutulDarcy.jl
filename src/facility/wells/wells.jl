

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

function setup_well(g, K, reservoir_cells::AbstractVector;
        simple_well = true,
        N = missing,
        neighborship = N,
        perforation_cells_well = missing,
        reference_depth = nothing,
        cell_centers = nothing,
        skin = 0.0,
        radius = 0.1,
        grouting_thickness = 0.0, # In addition to the radius
        casing_thickness = 0.0, # How much of the radius is casing
        Kh = missing,
        WI = missing,
        WIth = missing,
        volumes = missing,
        thermal_conductivity = missing,
        material_thermal_conductivity = 0.0,
        thermal_conductivity_casing = 20.0,
        thermal_conductivity_grout = 2.3,
        casing_heat_capacity = 420.0,
        casing_density = 8000.0,
        volume_multiplier = 1.0,
        friction = 1e-4, # Old version of kwarg for roughness
        roughness = friction,
        net_to_gross = missing,
        cell_radius = missing,
        well_cell_centers = missing,
        use_top_node = missing,
        dir = :z,
        drainage_radius = NaN,
        kwarg...
    )
    is_3d = dim(g) == 3
    n = length(reservoir_cells)
    # Make sure these are cell indices
    reservoir_cells = map(i -> cell_index(g, i), reservoir_cells)
    # Set up well itself
    if simple_well
        W = SimpleWell(reservoir_cells;
            kwarg...
        )
        perf_to_wellcell_index = [1]
    else
        if ismissing(neighborship)
            if ismissing(use_top_node)
                if is_3d && reference_depth isa Real
                    top_res_cell = reservoir_cells[1]
                    use_top_node = !(reference_depth ≈ cell_centers[3, top_res_cell])
                else
                    use_top_node = false
                end
            end
            W = MultiSegmentWell(reservoir_cells;
                top_node = use_top_node,
                kwarg...
            )
            perf_to_wellcell_index = collect(eachindex(reservoir_cells))
            if use_top_node
                pushfirst!(perf_to_wellcell_index, 1)
            end
        else
            # Well has actual topology
            !ismissing(perforation_cells_well) || error("Must provide perforation_cells_well if neighborship is provided.")
            W = MultiSegmentWell(neighborship, reservoir_cells, perforation_cells_well;
                kwarg...
            )
            perf_to_wellcell_index = perforation_cells_well
        end
    end
    treat_defaulted(x) = x
    treat_defaulted(::Missing) = NaN
    treat_defaulted(::Nothing) = NaN

    WI = treat_defaulted(WI)
    WIth = treat_defaulted(WIth)
    Kh = treat_defaulted(Kh)

    if isnothing(cell_centers)
        geometry = tpfv_geometry(g)
        cell_centers = geometry.cell_centroids
    end
    perforation_centers = cell_centers[:, reservoir_cells]
    if ismissing(well_cell_centers)
        if !simple_well && !ismissing(neighborship)
            error("Must provide well_cell_centers for multisegment wells when neighborship is provided.")
        end
        well_cell_centers = perforation_centers[:, perf_to_wellcell_index]
    else
        well_cell_centers = copy(well_cell_centers)
    end
    closest_reservoir_cell_to_well_cell = Int[]
    for i in axes(well_cell_centers, 2)
        dists = vec(norm.(eachcol(perforation_centers .- well_cell_centers[:, i]), 2))
        push!(closest_reservoir_cell_to_well_cell, findmin(dists)[2])
    end
    # Set reference depth if provided
    if reference_depth isa Real && is_3d
        well_cell_centers[3, 1] = reference_depth
    end

    Wdomain = DataDomain(W)
    c = Cells()
    p = Perforations()
    f = Faces()

    function cell_height(dir_i, cell)
        if dir_i isa Symbol
            Δ = cell_dims(g, cell)
            d_index = findfirst(isequal(dir_i), [:x, :y, :z])
            h = Δ[d_index]
        else
            h = norm(dir_i, 2)
        end
    end
    # ## Perforations

    Wdomain[:Kh, p] = Kh
    Wdomain[:skin, p] = skin
    Wdomain[:perforation_radius, p] = radius
    Wdomain[:well_index, p] = WI
    Wdomain[:perforation_centroids, p] = perforation_centers
    Wdomain[:drainage_radius, p] = drainage_radius
    if dir isa Symbol
        dir = fill(dir, n)
    end
    Wdomain[:perforation_direction, p] = dir

    direction_expanded = Wdomain[:perforation_direction, p]

    perf_height = map((i, cell) -> cell_height(direction_expanded[i], cell), 1:n, reservoir_cells)
    Wdomain[:cell_length, c] = perf_height[closest_reservoir_cell_to_well_cell]
    Wdomain[:cell_dims, p] = map(c -> peaceman_cell_dims(g, c), reservoir_cells)

    # Length of cell in direction of well - used for volumes of nodes
    if ismissing(cell_radius)
        cell_radius = Wdomain[:perforation_radius, p][1]
    end
    Wdomain[:radius, c] = cell_radius
    if !ismissing(volumes)
        length(volumes) == number_of_cells(W) || error("Must provide one volume per well cell ($(length(volumes)) provided, $(number_of_cells(W)) well cells).")
        Wdomain[:volume_override, c] = volumes
    end
    # Centers
    Wdomain[:cell_centroids, c] = well_cell_centers
    # Geometry
    Wdomain[:volume_multiplier, c] = volume_multiplier
    Wdomain[:casing_thickness, c] = casing_thickness
    Wdomain[:grouting_thickness, c] = grouting_thickness
    Wdomain[:thermal_conductivity_casing, c] = thermal_conductivity_casing
    Wdomain[:thermal_conductivity_grout, c] = thermal_conductivity_grout
    Wdomain[:casing_heat_capacity, c] = casing_heat_capacity
    Wdomain[:casing_density, c] = casing_density

    # ## Thermal well props
    # ### Perforations
    Wdomain[:thermal_well_index, p] = WIth
    # Perforation properties taken from reservoir
    perf_subset(x::AbstractVector) = x[reservoir_cells]
    perf_subset(x::AbstractMatrix) = x[:, reservoir_cells]

    if ismissing(net_to_gross)
        ntg = 1.0
    else
        ntg = perf_subset(net_to_gross)
    end
    Wdomain[:net_to_gross, p] = ntg
    Wdomain[:permeability, p] = perf_subset(K)
    if !ismissing(thermal_conductivity)
        Wdomain[:thermal_conductivity, p] = perf_subset(thermal_conductivity)
    end

    if !simple_well
        Wdomain[:roughness, f] = roughness
        Wdomain[:material_thermal_conductivity, f] = material_thermal_conductivity
    end

    return Wdomain
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
    if haskey(D, :net_to_gross)
        ntg = D[:net_to_gross]
    else
        ntg = missing
    end
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
    return setup_well(g, K, reservoir_cells;
        thermal_conductivity = Λ,
        cell_centers = cell_centers,
        net_to_gross = ntg,
        kwarg...
    )
end

function setup_well(g, K, reservoir_cell::Union{Int, Tuple, NamedTuple}; kwarg...)
    return setup_well(g, K, [reservoir_cell]; kwarg...)
end

function setup_well_from_trajectory(D::DataDomain, traj; traj_arg = NamedTuple(), kwarg...)
    G = D |> physical_representation |> UnstructuredMesh
    cells, extra = Jutul.find_enclosing_cells(G, traj; extra_out = true, traj_arg...)
    dir = Vector.(extra[:direction].*extra[:lengths])
    return setup_well(D, cells; dir = dir, kwarg...)
end

function map_well_nodes_to_reservoir_cells(w::MultiSegmentWell, reservoir::Union{DataDomain, Missing} = missing)
    # TODO: Try to more or less match it up cell by cell. Could be
    # improved...
    c = zeros(Int, number_of_cells(w))
    c[w.perforations.self] .= w.perforations.reservoir
    for i in eachindex(c)
        if i == firstindex(c)
            continue
        end
        if c[i] == 0
            c[i] = c[i-1]
        end
    end
    for i in reverse(eachindex(c))
        if i == lastindex(c)
            continue
        end
        if c[i] == 0
            c[i] = c[i+1]
        end
    end
    @assert all(x -> x > 0, c)
    return c
end

function map_well_nodes_to_reservoir_cells(w::DataDomain, reservoir::Union{DataDomain, Missing} = missing)
    return map_well_nodes_to_reservoir_cells(physical_representation(w), reservoir)
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

function domain_fluid_volume(d::DataDomain, grid::WellDomain)
    return domain_bulk_volume(d, grid, outer_boundary = :hole)
end

function domain_bulk_volume(d::DataDomain, grid::WellDomain; outer_boundary = :grouting)
    if haskey(d, :volume_override)
        vols = d[:volume_override, Cells()]
    else
        case_thickness = d[:casing_thickness, Cells()]
        grouting_thickness = d[:grouting_thickness, Cells()]
        mult = d[:volume_multiplier, Cells()]
        if grid isa MultiSegmentWell
            hole_radius = d[:radius, Cells()]
            r = well_bulk_volume_radius(hole_radius, case_thickness, grouting_thickness, outer_boundary = outer_boundary)
            L = d[:cell_length, Cells()]
            vols = mult.*(π .* r.^2 .* L)
        else
            # Simple wells are not segmented, so sum over perforations instead
            grid::SimpleWell
            ic = grid.perforations.self
            r = well_bulk_volume_radius(d[:perforation_radius, Perforations()], case_thickness[ic], grouting_thickness[ic], outer_boundary = outer_boundary)
            cdims = d[:cell_dims, Perforations()]
            dir = d[:perforation_direction, Perforations()]
            L = length_from_cell_dims.(cdims, dir)
            vols = only(mult)*sum(π .* r.^2 .* L)
        end
    end
    return vols
end

function well_bulk_volume_radius(r_h, r_c, r_g; outer_boundary::Symbol)
    if outer_boundary == :hole
        r = r_h
    elseif outer_boundary == :casing
        r = r_h .+ r_c
    elseif outer_boundary == :grouting
        r = r_h .+ r_c .+ r_g
    else
        error("Invalid outer_boundary: $outer_boundary, must be :hole, :casing or :grouting.")
    end
    return r
end

# Well segments
function get_neighborship(::SimpleWell)
    # No interior connections.
    return zeros(Int, 2, 0)
end

function number_of_cells(W::SimpleWell)
    return 1
end

function number_of_cells(W::MultiSegmentWell)
    return W.num_nodes
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
        ρ_w = state_well.PhaseMassDensities[ph, wc]
        s_w = state_well.Saturations[ph, wc]
        q_ph = s_w*ρ_w*ψ*λ_t
    else
        rc = conn.reservoir
        # Production
        if haskey(state_res, :PhaseMassMobilities)
            ρλ = state_res.PhaseMassMobilities[ph, rc]
        else
            ρ = state_res.PhaseMassDensities[ph, rc]
            λ = state_res.PhaseMobilities[ph, rc]
            ρλ = ρ*λ
        end
        q_ph = ρλ*ψ
    end
    return q_ph
end

Base.@propagate_inbounds function simple_well_perforation_flux!(out, sys::Union{ImmiscibleSystem, SinglePhaseSystem}, state_res, state_well, rhoS, conn)
    rc = conn.reservoir
    # Reservoir quantities
    ρ = state_res.PhaseMassDensities
    mob = state_res.PhaseMobilities
    nph = size(ρ, 1)
    ρλ_t = zero(eltype(mob))
    for ph in 1:nph
        ρλ_t += ρ[ph, rc]*mob[ph, rc]
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
            λ = mob[ph, rc]
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

function compute_well_cell_volumes(; radius = missing, lengths = missing, volumes = missing)
    if ismissing(volumes)
        !ismissing(radius) || error("Either radius or volumes must be provided")
        !ismissing(lengths) || error("Either lengths or volumes must be provided")
        n = length(lengths)
        T = promote_type(eltype(radius), eltype(lengths))
        volumes = zeros(T, n)
        for i in 1:n
            volumes[i] = π*radius^2*lengths[i]
        end
    end
    return volumes
end
