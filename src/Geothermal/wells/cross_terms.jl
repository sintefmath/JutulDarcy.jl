## Abstract types for BTES cross terms

abstract type AbstractBTESCT <: Jutul.AdditiveCrossTerm end

Jutul.symmetry(::AbstractBTESCT) = Jutul.CTSkewSymmetry()

function Jutul.cross_term_entities(ct::AbstractBTESCT, eq::ConservationLaw, model)
    return ct.btes_cells
end

function Jutul.cross_term_entities_source(ct::AbstractBTESCT, eq::ConservationLaw, model)
    return ct.btes_cells
end

## Cross terms types for BTES wells (supply-return mass/energy, grout-grout energy)

struct BTESWellSupplyToReturnMassCT <: AbstractBTESCT
    btes_cells::Vector{Int64}
end

struct BTESWellSupplyToReturnEnergyCT <: AbstractBTESCT
    btes_cells::Vector{Int64}
end

struct BTESWellGroutEnergyCT <: AbstractBTESCT
    WIth_grout::Vector{Float64}
    btes_cells::Vector{Int64}
end

## Cross term calculations for BTES wells

function JutulDarcy.update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::BTESWellSupplyToReturnMassCT, eq, dt, ldisc = JutulDarcy.local_discretization(ct, i))

    # Unpack properties
    sys = JutulDarcy.flow_system(model_t.system)
    @inbounds begin 
        btes_cell = ct.btes_cells[i]
        nph = number_of_phases(sys)
        for ph in 1:nph
            q_ph = btes_supply_return_massflux(state_s, state_t, btes_cell, ph)
            out[ph] = q_ph
        end

    end

end

function JutulDarcy.update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::BTESWellSupplyToReturnEnergyCT, eq, dt, ldisc = JutulDarcy.local_discretization(ct, i))

    # Unpack properties
    sys = JutulDarcy.flow_system(model_t.system)
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

function JutulDarcy.update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::BTESWellGroutEnergyCT, eq, dt, ldisc = JutulDarcy.local_discretization(ct, i))

    # Unpack properties
    sys = JutulDarcy.flow_system(model_t.system)
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
