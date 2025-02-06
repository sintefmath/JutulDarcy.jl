function setup_btes_well(D::DataDomain, reservoir_cells;
    btes_type = :simple,
    name = :BTES, kwarg...)

    if btes_type == :simple
        return setup_btes_well_simple(D, reservoir_cells; name = name, kwarg...)

    elseif btes_type == :u1
        return setup_btes_well_u1(D, reservoir_cells; name = name, kwarg...)

    else
        error("Unknown BTES type: $btes_type")
    end

end

function setup_vertical_btes_well(D::DataDomain, i, j; heel = 1, toe = missing, kwarg...)

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
    setup_btes_well(D, reservoir_cells; kwarg...)

end

function setup_btes_well_simple(D::DataDomain, reservoir_cells;
    name = :BTES,
    WIth_grout = missing,
    kwarg...)

    if ismissing(WIth_grout)
        WIth_grout = zeros(length(reservoir_cells))
    end
    @assert length(WIth_grout) == length(reservoir_cells)

    args = (
        WI = 0.0,
        extra_perforation_props = (WIth_grout = WIth_grout, ),
        simple_well = false,
        type = :btes
    )

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
    radius_grout = 60e-3,
    radius_pipe_inner = 20e-3,
    radius_pipe_outer = radius_pipe_inner + 5e-3,
    pipe_spacing = radius_grout,
    thermal_conductivity_grout = 2.0,
    thermal_conductivity_pipe = 2.0,
    friction = 1e-4,
    dir = :z,
    kwarg...)

    # Set up connectivity
    #--------------------------------------------------------------------------#
    nc_pipe = length(reservoir_cells)
    nc_grout = length(reservoir_cells)

    pipe_cells = (1:nc_pipe)
    grout_cells = (1:nc_grout) .+ nc_pipe

    pipe_to_pipe = vcat(pipe_cells[1:end-1]', pipe_cells[2:end]')
    pipe_to_grout = vcat(pipe_cells', grout_cells')
    grout_to_grout = vcat(grout_cells[1:end-1]', grout_cells[2:end]')

    N = hcat(pipe_to_pipe, pipe_to_grout)
    # N = hcat(pipe_to_pipe, pipe_to_grout, grout_to_grout)
    #--------------------------------------------------------------------------#

    # Compute BTES segment properties
    #--------------------------------------------------------------------------#
    function get_entry(x::AbstractVector, i)
        return x[i]
    end
    function get_entry(x, i)
        return x
    end
    mesh = physical_representation(D)

    btes_type = BTESTypeU1()

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
        vol_p[i], vol_g[i], L = btes_volume(btes_type, mesh, c, d, rg, rp_in, rp_out)
        λpg[i], λgr[i], λgg[i] = btes_thermal_conductivity(btes_type, rg, rp_in, rp_out, w, L, λg, λp)
    end
    volumes = vcat(vol_p, vol_g)
    #--------------------------------------------------------------------------#

    perforation_cells = collect(grout_cells)

    # Set up segment flow models
    #--------------------------------------------------------------------------#
    centers = repeat(cell_centers[:, reservoir_cells], 1, 2)
    reference_depth = centers[3, 1]
    accumulator_center = [centers[1, 1], centers[2, 1], reference_depth]
    ext_centers = hcat(accumulator_center, centers)

    segment_models = Vector{Any}()
    N0 = copy(N)
    N = hcat([1,2], N.+1)
    # println(N)

    # println(λgr)
    nseg = size(N,2)
    for seg in 1:nseg
        l, r = N[:, seg]
        L = norm(ext_centers[:, l] - ext_centers[:, r], 2)
        if seg <= nc_pipe+1
            Do, Di = 2*radius_pipe_inner, 0.0
        else
            Do, Di = 2*radius_pipe_inner, 2*radius_pipe_inner*0.99
        end
        Δp = SegmentWellBoreFrictionHB(L, friction, Do; D_inner = Di)
        push!(segment_models, Δp)
    end

    dz = cell_centers[3, reservoir_cells] .- reference_depth
    
    nr = length(reservoir_cells)

    material_thermal_conductivity = zeros(nseg)
    material_thermal_conductivity[grout_cells] .= λpg

    WI = fill(0.0, nr)
    supply_well = MultiSegmentWell(reservoir_cells, volumes, centers;
        type = :btes,
        N = N0,
        WI = WI,
        WIth = λgr,
        extra_perforation_props = (WIth_grout = λgg, ),
        material_thermal_conductivity = material_thermal_conductivity,
        dz = dz,
        name = Symbol(name, "_supply"),
        perforation_cells = perforation_cells,
        segment_models = segment_models,
    )

    return_well = MultiSegmentWell(reservoir_cells, volumes, centers;
        type = :btes,
        N = N0,
        WI = WI,
        WIth = λgr,
        extra_perforation_props = (WIth_grout = λgg, ),
        material_thermal_conductivity = material_thermal_conductivity,
        dz = dz,
        name = Symbol(name, "_return"),
        perforation_cells = perforation_cells,
        segment_models = segment_models,
    )


    return supply_well, return_well

end

abstract type AbstractBTESType end

struct BTESTypeU1 <: AbstractBTESType end

function btes_volume(type::BTESTypeU1, g, reservoir_cell, dir, radius_grout, radius_pipe_inner, radius_pipe_outer)

    # Get BTES segment length
    
    Δ = cell_dims(g, reservoir_cell)
    dir_index = findfirst(isequal(dir), [:x, :y, :z])
    L = Δ[dir_index]

    # Compute pipe and grout volume
    rg, rp_in, rp_out = radius_grout, radius_pipe_inner, radius_pipe_outer
    vol_p = π*(rp_out^2 - rp_in^2)*L
    vol_g = π*rg^2*L/2 - vol_p

    return vol_p, vol_g, L

end

function btes_thermal_conductivity(type::BTESTypeU1, 
    radius_grout, radius_pipe_inner, radius_pipe_outer, pipe_spacing, length,
    thermal_conductivity_pipe, thermal_conductivity_grout)

    # Conveient short-hand notation
    #--------------------------------------------------------------------------#
    rg, rp_in, rp_out, w, L = 
        radius_grout, radius_pipe_inner, radius_pipe_outer, pipe_spacing, length
    λg, λp = thermal_conductivity_grout, thermal_conductivity_pipe
    #--------------------------------------------------------------------------#

    # Compute thermal resistances
    #--------------------------------------------------------------------------#
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

    Rar = acosh((2*w^2 - dp_out^2)/dp_out^2)/(2*π*λg)
    Rgg = 2*Rgr*(Rar - 2*x*Rg)/(2*Rgr - Rar + 2*x*Rg)
    #--------------------------------------------------------------------------#

    # Compute thermal transmissibilities
    #--------------------------------------------------------------------------#
    λpg = L*2*π*rp_in/Rpg
    λgr = L*π*rg/Rgr
    λgg = L*2*rg/Rg
    #--------------------------------------------------------------------------#

    return λpg, λgr, λgg

end

struct MaterialThermalConductivities <: ScalarVariable end
Jutul.variable_scale(::MaterialThermalConductivities) = 1e-10
Jutul.minimum_value(::MaterialThermalConductivities) = 0.0
Jutul.associated_entity(::MaterialThermalConductivities) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::MaterialThermalConductivities, symb)
    if haskey(data_domain, :material_thermal_conductivity, Faces())
        T = copy(data_domain[:material_thermal_conductivity])
    else
        error(":material_thermal_conductivity or :material_thermal_conductivity symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

struct SegmentWellBoreFrictionDarcy{R}
    L::R
    area::R
    perm::R
    function SegmentWellBoreFrictionDarcy(L, perm)
        new{typeof(L)}(L, area, perm)
    end
end

struct BTESWellSupplyToReturnMassCT <: Jutul.AdditiveCrossTerm
    btes_cells::Vector{Int64}
end

struct BTESWellSupplyToReturnEnergyCT <: Jutul.AdditiveCrossTerm
    btes_cells::Vector{Int64}
end

struct BTESWellGroutEnergyCT <: Jutul.AdditiveCrossTerm
    WIth_grout::Vector{Float64}
    btes_cells::Vector{Int64}
end

function Jutul.cross_term_entities(ct::BTESWellSupplyToReturnMassCT, eq::ConservationLaw, model)
    return ct.btes_cells
end

function Jutul.cross_term_entities_source(ct::BTESWellSupplyToReturnMassCT, eq::ConservationLaw, model)
    return ct.btes_cells
end

function Jutul.cross_term_entities(ct::BTESWellSupplyToReturnEnergyCT, eq::ConservationLaw, model)
    return ct.btes_cells
end

function Jutul.cross_term_entities_source(ct::BTESWellSupplyToReturnEnergyCT, eq::ConservationLaw, model)
    return ct.btes_cells
end

function Jutul.cross_term_entities(ct::BTESWellGroutEnergyCT, eq::ConservationLaw, model)
    return ct.btes_cells
end

function Jutul.cross_term_entities_source(ct::BTESWellGroutEnergyCT, eq::ConservationLaw, model)
    return ct.btes_cells
end

Jutul.symmetry(::BTESWellSupplyToReturnMassCT) = Jutul.CTSkewSymmetry()
Jutul.symmetry(::BTESWellSupplyToReturnEnergyCT) = Jutul.CTSkewSymmetry()
Jutul.symmetry(::BTESWellGroutEnergyCT) = Jutul.CTSkewSymmetry()

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::BTESWellSupplyToReturnMassCT, eq, dt, ldisc = local_discretization(ct, i))

    # Unpack properties
    sys = flow_system(model_t.system)
    @inbounds begin 
        btes_cell = ct.btes_cells[i]
        nph = number_of_phases(sys)
        for ph in 1:nph
            q_ph = btes_supply_return_massflux(state_s, state_t, btes_cell, ph)
            out[ph] = q_ph
        end

    end

end

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::BTESWellSupplyToReturnEnergyCT, eq, dt, ldisc = local_discretization(ct, i))

    # Unpack properties
    sys = flow_system(model_t.system)
    @inbounds begin 
        btes_cell = ct.btes_cells[i]
        nph = number_of_phases(sys)
        q = 0.0
        for ph in 1:nph
            q_ph = btes_supply_return_massflux(state_s, state_t, btes_cell, ph)
            h_ph = state_s.FluidEnthalpy[ph,btes_cell]
            q += q_ph.*h_ph
        end
    end
    out[] = q

end

Base.@propagate_inbounds function btes_supply_return_massflux(state_supply, state_return, cell, ph)

    p_s = state_supply.Pressure[cell]
    p_t = state_return.Pressure[cell]
    dp = p_s - p_t
    
    T = 1.0e-10
    Ψ = -T.*dp
    
    ρ = state_supply.PhaseMassDensities[ph,cell]
    s = state_supply.Saturations[ph,cell]
    μ = state_supply.PhaseViscosities[ph,cell]

    q_ph = s.*ρ./μ.*Ψ

    return q_ph
        
end

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::BTESWellGroutEnergyCT, eq, dt, ldisc = local_discretization(ct, i))

    # Unpack properties
    sys = flow_system(model_t.system)
    @inbounds begin 
        btes_cell = ct.btes_cells[i]
        λ = ct.WIth_grout[i]
        T_s = state_s.Temperature[btes_cell]
        T_t = state_t.Temperature[btes_cell]
        λ = ct.WIth_grout[i]
        q = -λ.*(T_t - T_s)
    end

    out[] = q

end