export ReservoirFromFractureFlowCT

abstract type AbstractReservoirFromFractureCT <: Jutul.AdditiveCrossTerm end

struct ReservoirFromFractureFlowCT{I<:AbstractVector, T<:AbstractVector} <: AbstractReservoirFromFractureCT
    reservoir_cells::I
    fracture_cells::I
    connection_strength::T
end

function ReservoirFromFractureFlowCT(reservoir_cells::Vector{Int}, fracture_cells::Vector{Int}, connection_strength::Vector{Float64})
    return ReservoirFromFractureFlowCT{Vector{Int}, Vector{Float64}}(reservoir_cells, fracture_cells, connection_strength)
end

function Base.show(io::IO, d::ReservoirFromFractureFlowCT)
    n = length(d.reservoir_cells)
    print(io, "ReservoirFromFractureFlowCT ($n connections)")
end

Jutul.symmetry(::AbstractReservoirFromFractureCT) = Jutul.CTSkewSymmetry()

function Jutul.cross_term_entities(ct::AbstractReservoirFromFractureCT, eq::ConservationLaw, model)
    return ct.reservoir_cells
end

function Jutul.cross_term_entities_source(ct::AbstractReservoirFromFractureCT, eq::ConservationLaw, model)
    return ct.fracture_cells
end

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::ReservoirFromFractureFlowCT, eq, dt, ldisc = local_discretization(ct, i))
    
    # Target (Matrix)
    rc = ct.reservoir_cells[i]
    # Source (Fracture)
    fc = ct.fracture_cells[i]
    # Transmissibility
    T_rf = ct.connection_strength[i]

    p_m = state_t.Pressure[rc]
    p_f = state_s.Pressure[fc]
    # Driving force: Fracture pressure - Matrix pressure
    Δp = p_f - p_m

    # Properties
    rho_m = state_t.PhaseMassDensities
    rho_f = state_s.PhaseMassDensities
    
    mob_m = state_t.PhaseMobilities
    mob_f = state_s.PhaseMobilities

    nph = size(out, 1)
    @inbounds for ph in 1:nph
        # Upwinding
        if Δp > 0
            # Fracture -> Matrix
            ρ = rho_f[ph, fc]
            λ = mob_f[ph, fc]
        else
            # Matrix -> Fracture
            ρ = rho_m[ph, rc]
            λ = mob_m[ph, rc]
        end
        out[ph] = T_rf * ρ * λ * Δp
    end
end

function Jutul.subcrossterm(ct::ReservoirFromFractureFlowCT, ctp, m_t, m_s, map_res::Jutul.FiniteVolumeGlobalMap, map_frac::Jutul.FiniteVolumeGlobalMap, partition)
    rc = ct.reservoir_cells[ctp]
    fc = ct.fracture_cells[ctp]
    T = ct.connection_strength[ctp]

    rc_local = map(c -> Jutul.local_cell(c, map_res), rc)
    fc_local = map(c -> Jutul.local_cell(c, map_frac), fc)

    return ReservoirFromFractureFlowCT(rc_local, fc_local, T)
end