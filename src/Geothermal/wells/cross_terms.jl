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
    sys = model_t.system
    @inbounds begin 
        supply_node = ct.supply_nodes[i]
        return_node = ct.return_nodes[i]
        q = btes_supply_return_massflux(sys, state_s, state_t, supply_node, return_node)
        for j in eachindex(out)
            out[j] = q[j]
        end

    end

end

function JutulDarcy.update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::ClosedLoopSupplyToReturnEnergyCT, eq, dt, ldisc = JutulDarcy.local_discretization(ct, i))

    # Unpack properties
    sys = model_t.system
    @inbounds begin 
        supply_node = ct.supply_nodes[i]
        return_node = ct.return_nodes[i]
        q = btes_supply_return_heatflux(sys, state_s, state_t, supply_node, return_node)
    end
    out[] = q

end

Base.@propagate_inbounds function btes_supply_return_massflux(system::JutulSystem, state_supply, state_return, supply_node, return_node)

    T = 1.0e-10
    p_s = state_supply.Pressure[supply_node]
    p_t = state_return.Pressure[return_node]
    dp = p_s - p_t
    Ψ = -T*dp

    nph = number_of_phases(system)
    S = Base.promote_type(typeof(Ψ), typeof(p_s), typeof(p_t))
    q = zeros(S, nph)
    if dp > 0
        state, node = state_supply, supply_node
    else
        state, node = state_return, return_node
    end
    for ph in 1:nph
        ρ = state.PhaseMassDensities[ph, node]
        s = state.Saturations[ph, node]
        μ = state.PhaseViscosities[ph, node]
        q[ph] = (s*ρ/μ)*Ψ
    end

    return q
end

Base.@propagate_inbounds function btes_supply_return_massflux(system::JutulDarcy.AbstractCompositionalSystemLV, state_supply, state_return, supply_node, return_node)

    T = 1.0e-10
    p_s = state_supply.Pressure[supply_node]
    p_t = state_return.Pressure[return_node]
    dp = p_s - p_t
    Ψ = -T*dp

    nph = number_of_phases(system)
    ncomp = JutulDarcy.number_of_components(system)
    S = Base.promote_type(typeof(Ψ), typeof(p_s), typeof(p_t))
    q = zeros(S, ncomp)
    if dp > 0
        state, node = state_supply, supply_node
    else
        state, node = state_return, return_node
    end
    λ_t = 0.0
    rho_mix = 0.0
    for ph in 1:nph
        λ_t += state.Saturations[ph, node]/state.PhaseViscosities[ph, node]
        rho_mix += state.Saturations[ph, node]*state.PhaseMassDensities[ph, node]
    end
    mass_t = 0.0
    for i = 1:ncomp
        mass_t += state.TotalMasses[i, node]
    end
    for i = 1:ncomp
        f = state.TotalMasses[i, node]/mass_t
        q[i] = f*λ_t*rho_mix*Ψ
    end

    return q
end

Base.@propagate_inbounds function btes_supply_return_heatflux(system::JutulSystem, state_supply, state_return, supply_node, return_node)
    q_ph = btes_supply_return_massflux(system, state_supply, state_return, supply_node, return_node)
    nph = number_of_phases(system)
    q = 0.0
    for ph in 1:nph
        if q_ph[ph] < 0
            h_ph = state_supply.FluidEnthalpy[ph, supply_node]
        else
            h_ph = state_return.FluidEnthalpy[ph, return_node]
        end
        q += q_ph[ph]*h_ph
    end
   
    return q
end

Base.@propagate_inbounds function btes_supply_return_heatflux(system::JutulDarcy.AbstractCompositionalSystemLV, state_supply, state_return, supply_node, return_node)
    q_comp = btes_supply_return_massflux(system, state_supply, state_return, supply_node, return_node)
    ncomp = JutulDarcy.number_of_components(system)
    q = 0.0
    has_enthalpy = haskey(state_supply, :Enthalpy) && haskey(state_return, :Enthalpy)
    function enthalpy(comp, state, node)
        if has_enthalpy
            return state.Enthalpy[node]
        else
            h_l = state.FluidEnthalpy[1, node]
            h_v = state.FluidEnthalpy[2, node]
            S_l = state.Saturations[1, node]
            S_v = state.Saturations[2, node]
            return S_l*h_l + S_v*h_v
        end
    end
    for comp in 1:ncomp
        if q_comp[comp] < 0
            h_comp = enthalpy(comp, state_supply, supply_node)
        else
            h_comp = enthalpy(comp, state_return, return_node)
        end
        q += q_comp[comp]*h_comp
    end
   
    return q
end

function JutulDarcy.update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::BTESWellGroutEnergyCT, eq, dt, ldisc = JutulDarcy.local_discretization(ct, i))

    # Unpack properties
    sys = model_t.system
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
