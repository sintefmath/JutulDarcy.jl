using LinearAlgebra

## Utility functions for setting up BTES wells
function setup_btes_well(D::DataDomain, reservoir_cells;
    btes_type = :simple,
    name = :BTES, kwarg...)

    if btes_type == :simple
        # Simple BTES well with direct pipe/reservoir thermal communication
        return setup_btes_well_simple(D, reservoir_cells; name = name, kwarg...)

    elseif btes_type == :u1
        # BTES well with U1 type geometry, accounting for accumulation of
        # thermal energy in grouting
        return setup_btes_well_u1(D, reservoir_cells; name = name, kwarg...)

    else
        # Unknown BTES type
        error("Unknown BTES type: $btes_type")
    end

end

function setup_vertical_btes_well(D::DataDomain, i, j;
    heel = 1, toe = missing, kwarg...)

    # Get reservoir cells from ijk-indices
    g = physical_representation(D)
    if ismissing(toe)
        toe = grid_dims_ijk(g)[3]
    end
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
    # Set up BTES well
    setup_btes_well(D, reservoir_cells; kwarg...)

end

function setup_btes_well_simple(D::DataDomain, reservoir_cells;
    name = :BTES,
    WIth_grout = missing,
    radius_pipe = 16e-3 - 2.9e-3,
    kwarg...)

    # Set thermal conductivities between supply and return
    if ismissing(WIth_grout)
        WIth_grout = zeros(length(reservoir_cells))
    end
    @assert length(WIth_grout) == length(reservoir_cells)

    # Common properties
    args = (
        WI = 0.0,
        extra_perforation_props = (WIth_grout = WIth_grout, ),
        radius = radius_pipe,
        end_nodes = [length(reservoir_cells)+1],
        simple_well = false,
        type = :closed_loop
    )

    # Set up supply and return wells
    supply_well = setup_well(D::DataDomain, reservoir_cells;
        name = Symbol(name, "_supply"),
        args..., kwarg...)
    return_well = setup_well(D::DataDomain, reservoir_cells;
        name = Symbol(name, "_return"),
        args..., kwarg...)

    return supply_well, return_well

end

function setup_btes_well_u1(D::DataDomain, reservoir_cells;
    cell_centers = D[:cell_centroids],
    name = :BTES,
    radius_grout = 65e-3,
    radius_pipe_outer = 16e-3,
    radius_pipe_inner = radius_pipe_outer - 2.9e-3,
    pipe_spacing = 60e-3,
    thermal_conductivity_grout = 2.3,
    thermal_conductivity_pipe = 0.38,
    # thermal_conductivity_pipe = 2.0,
    heat_capacity_grout = 420.0,
    density_grout = 1500.0,
    friction = 1e-4,
    dir = :z,
    kwarg...)

    ## Set up connectivity
    nc_pipe = length(reservoir_cells)
    nc_grout = length(reservoir_cells)

    pipe_cells = (1:nc_pipe)
    grout_cells = (1:nc_grout) .+ nc_pipe

    pipe_to_pipe = vcat(pipe_cells[1:end-1]', pipe_cells[2:end]')
    pipe_to_grout = vcat(pipe_cells', grout_cells')
    grout_to_grout = vcat(grout_cells[1:end-1]', grout_cells[2:end]')

    N = hcat(pipe_to_pipe, pipe_to_grout)

    ## Compute BTES segment properties
    btes_type = BTESTypeU1()
    mesh = physical_representation(D)

    function get_entry(x::AbstractVector, i)
        return x[i]
    end
    function get_entry(x, i)
        return x
    end

    vol_p = zeros(nc_pipe)
    vol_g = zeros(nc_pipe)
    λpg = zeros(nc_pipe)
    λgr = zeros(nc_grout)
    λgg = zeros(nc_grout)

    λg, λp = thermal_conductivity_grout, thermal_conductivity_pipe

    for (i, c) = enumerate(reservoir_cells)
        # Get segment properties
        rg = get_entry(radius_grout, i)
        rp_in = get_entry(radius_pipe_inner, i)
        rp_out = get_entry(radius_pipe_outer, i)
        w = get_entry(pipe_spacing, i)
        d = (dir isa Symbol) ? dir : dir[i]
        # Compute thermal conductivities and volumes
        vol_p[i], vol_g[i], L = btes_volume(
            btes_type, mesh, c, d, rg, rp_in, rp_out)
        λpg[i], λgr[i], λgg[i] = btes_thermal_conductivity(
            btes_type, rg, rp_in, rp_out, w, L, λg, λp)
    end
    volumes = vcat(vol_p, vol_g)

    ## Set up segment flow models
    segment_models = Vector{Any}()

    # Set centers and depths
    centers = repeat(cell_centers[:, reservoir_cells], 1, 2)
    reference_depth = centers[3, 1]
    accumulator_center = [centers[1, 1], centers[2, 1], reference_depth]
    ext_centers = hcat(accumulator_center, centers)

    # Add top node
    N0 = copy(N)
    N = hcat([1,2], N.+1)
    nseg = size(N,2)

    # Set material thermal conducivities
    nr = length(reservoir_cells)

    for seg in 1:nseg
        l, r = N[:, seg]
        L = norm(ext_centers[:, l] - ext_centers[:, r], 2)
        if seg <= nc_pipe
            # Pipe segments use standard wellbore friction model
            Do, Di = 2*radius_pipe_inner, 0.0
            seg_model = SegmentWellBoreFrictionHB(L, friction, Do; D_inner = Di)
        else
            seg_model = JutulDarcy.ClosedSegment()
        end
        push!(segment_models, seg_model)
    end
    dz = cell_centers[3, reservoir_cells] .- reference_depth

    material_thermal_conductivity = zeros(nseg)
    material_thermal_conductivity[grout_cells] .= λpg

    nc = nc_pipe + nc_grout + 1
    void_fraction = ones(nc)
    void_fraction[grout_cells .+ 1] .= 0.0

    ## Set up supply and return wells
    args = (
        type = :closed_loop,
        N = N0,
        WI = fill(0.0, nr),
        WIth = λgr,
        extra_perforation_props = (WIth_grout = λgg, ),
        material_thermal_conductivity = material_thermal_conductivity,
        material_heat_capacity = fill(heat_capacity_grout, nc),
        material_density = fill(density_grout, nc),
        void_fraction = void_fraction,
        dz = dz,
        perforation_cells = collect(grout_cells.+1),
        end_nodes = [nc_pipe+1],
        segment_models = segment_models,
    )

    supply_well = MultiSegmentWell(reservoir_cells, volumes, centers;
        name = Symbol(name, "_supply"), args...)
    return_well = MultiSegmentWell(reservoir_cells, volumes, centers;
        name = Symbol(name, "_return"), args...)

    return supply_well, return_well

end

## Utility functions

# Conveience types for multiple dispatching
abstract type AbstractBTESType end
struct BTESTypeU1 <: AbstractBTESType end

function btes_volume(type::BTESTypeU1, g, reservoir_cell, dir, radius_grout, radius_pipe_inner, radius_pipe_outer)

    # Get BTES segment length
    Δ = cell_dims(g, reservoir_cell)
    dir_index = findfirst(isequal(dir), [:x, :y, :z])
    L = Δ[dir_index]

    # Compute pipe and grout volume
    rg, rp_in, rp_out = radius_grout, radius_pipe_inner, radius_pipe_outer
    vol_p = π*rp_in^2*L
    vol_g = π*rg^2*L/2 - vol_p

    return vol_p, vol_g, L

end

function btes_thermal_conductivity(type::BTESTypeU1, 
    radius_grout, radius_pipe_inner, radius_pipe_outer, pipe_spacing, length,
    thermal_conductivity_grout, thermal_conductivity_pipe)
    # Conveient short-hand notation
    rg, rp_in, rp_out, w, L = 
        radius_grout, radius_pipe_inner, radius_pipe_outer, pipe_spacing, length
    λg, λp = thermal_conductivity_grout, thermal_conductivity_pipe

    ## Compute thermal resistances
    # Advection-dependent pipe thermal resistance
    Ra = 0.0 # TODO: Implement this
    # Conduction-dependent pipe thermal resistance
    Rc_a = log(rp_out/rp_in)/(2*π*λp)
    dg, dp_out = 2*rg, 2*rp_out
    x = log(sqrt(dg^2 + 2*dp_out^2)/(2*dp_out))/log(dg/(sqrt(2)*dp_out))
    # Grout thermal resistance
    Rg = acosh((dg^2 + dp_out^2 - w^2)/(2*dg*dp_out))/(2*π*λg)*(1.601 - 0.888*w/dg)
    # Conduction-dependent grout thermal resistance
    Rc_b = x*Rg
    # Combined thermal resistance of pipe and grout
    Rpg = Ra + Rc_a + Rc_b
    Rgr = (1-x)*Rg
    # Grout-to-grout thermal resistance
    Rar = acosh((2*w^2 - dp_out^2)/dp_out^2)/(2*π*λg)
    Rgg = 2*Rgr*(Rar - 2*x*Rg)/(2*Rgr - Rar + 2*x*Rg)

    # Compute thermal conducivities
    λpg = L*2*π*rp_in/Rpg
    λgr = L*π*rg/Rgr
    λgg = L*2*rg/Rgg

    return λpg, λgr, λgg

end
