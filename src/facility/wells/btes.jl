function setup_btes_well(D::DataDomain, reservoir_cells; name = :BTES, kwarg...)

    supply_well = setup_well(D::DataDomain, reservoir_cells;
        name = Symbol(name, "_supply"), WI = 0.0, simple_well = false, type = :btes, kwarg...)

    return_well = setup_well(D::DataDomain, reservoir_cells;
        name = Symbol(name, "_return"), WI = 0.0, simple_well = false, type = :btes, kwarg...)
   
    return supply_well, return_well

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

struct BTESWellSupplyToReturnMassCT <: Jutul.AdditiveCrossTerm
    btes_cells::Vector{Int64}
end

struct BTESWellSupplyToReturnEnergyCT <: Jutul.AdditiveCrossTerm
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

Jutul.symmetry(::BTESWellSupplyToReturnMassCT) = Jutul.CTSkewSymmetry()
Jutul.symmetry(::BTESWellSupplyToReturnEnergyCT) = Jutul.CTSkewSymmetry()

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


function btes_supply_return_massflux(state_supply, state_return, cell, ph)

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