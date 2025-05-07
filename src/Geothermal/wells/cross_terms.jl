## Abstract types for BTES cross terms

abstract type AbstractClosedLoopCT <: Jutul.AdditiveCrossTerm end

Jutul.symmetry(::AbstractClosedLoopCT) = Jutul.CTSkewSymmetry()

function Jutul.cross_term_entities(ct::AbstractClosedLoopCT, eq::ConservationLaw, model)
    return ct.return_nodes
end

function Jutul.cross_term_entities_source(ct::AbstractClosedLoopCT, eq::ConservationLaw, model)
    return ct.supply_nodes
end

## Cross terms types for BTES wells (supply-return mass/energy, grout-grout energy)

struct ClosedLoopSupplyToReturnMassCT <: AbstractClosedLoopCT
    supply_nodes::Vector{Int64}
    return_nodes::Vector{Int64}
    function ClosedLoopSupplyToReturnMassCT(supply_nodes::Vector{Int64}, return_nodes::Vector{Int64})
        if length(supply_nodes) != length(return_nodes)
            throw(ArgumentError("Supply and return nodes must have the same length"))
        end
        new(supply_nodes, return_nodes)
    end
end

struct ClosedLoopSupplyToReturnEnergyCT <: AbstractClosedLoopCT
    supply_nodes::Vector{Int64}
    return_nodes::Vector{Int64}
    function ClosedLoopSupplyToReturnEnergyCT(supply_nodes::Vector{Int64}, return_nodes::Vector{Int64})
        if length(supply_nodes) != length(return_nodes)
            throw(ArgumentError("Supply and return nodes must have the same length"))
        end
        new(supply_nodes, return_nodes)
    end
end

struct BTESWellGroutEnergyCT <: AbstractClosedLoopCT
    WIth_grout::Vector{Float64}
    supply_nodes::Vector{Int64}
    return_nodes::Vector{Int64}
    function BTESWellGroutEnergyCT(WIth_grout::Vector{Float64}, supply_nodes::Vector{Int64}, return_nodes::Vector{Int64})
        if length(supply_nodes) != length(return_nodes)
            throw(ArgumentError("Supply and return nodes must have the same length"))
        end
        new(WIth_grout, supply_nodes, return_nodes)
    end
end

## Cross term calculations for BTES wells

function JutulDarcy.update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::ClosedLoopSupplyToReturnMassCT, eq, dt, ldisc = JutulDarcy.local_discretization(ct, i))

    # Unpack properties
    sys = JutulDarcy.flow_system(model_t.system)
    @inbounds begin 
        supply_node = ct.supply_nodes[i]
        return_node = ct.return_nodes[i]
        nph = number_of_phases(sys)
        for ph in 1:nph
            q_ph = btes_supply_return_massflux(state_s, state_t, supply_node, return_node, ph)
            out[ph] = q_ph
        end

    end

end

function JutulDarcy.update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::ClosedLoopSupplyToReturnEnergyCT, eq, dt, ldisc = JutulDarcy.local_discretization(ct, i))

    # Unpack properties
    sys = JutulDarcy.flow_system(model_t.system)
    @inbounds begin 
        supply_node = ct.supply_nodes[i]
        return_node = ct.return_nodes[i]
        nph = number_of_phases(sys)
        q = 0.0
        for ph in 1:nph
            q_ph = btes_supply_return_massflux(state_s, state_t, supply_node, return_node, ph)
            h_ph = state_s.FluidEnthalpy[ph, supply_node]
            q += q_ph.*h_ph
        end
    end
    out[] = q

end

Base.@propagate_inbounds function btes_supply_return_massflux(state_supply, state_return, supply_node, return_node, ph)

    p_s = state_supply.Pressure[supply_node]
    p_t = state_return.Pressure[return_node]
    dp = p_s - p_t

    T = 1.0e-10
    Ψ = -T.*dp

    ρ = state_supply.PhaseMassDensities[ph, supply_node]
    s = state_supply.Saturations[ph, supply_node]
    μ = state_supply.PhaseViscosities[ph, supply_node]

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
        supply_node = ct.supply_nodes[i]
        return_node = ct.return_nodes[i]
        λ = ct.WIth_grout[i]
        T_s = state_s.Temperature[supply_node]
        T_t = state_t.Temperature[return_node]
        λ = ct.WIth_grout[i]
        q = -λ.*(T_t - T_s)
    end

    out[] = q
end
