abstract type AbstractMatrixFromFractureCT <: Jutul.AdditiveCrossTerm end

struct MatrixFromFractureFlowCT{I<:AbstractVector, T<:AbstractVector} <: AbstractMatrixFromFractureCT
    matrix_cells::I
    fracture_cells::I
    connection_strength::T
end

function MatrixFromFractureFlowCT(matrix_cells::Vector{Int}, fracture_cells::Vector{Int}, connection_strength::Vector{Float64})
    return MatrixFromFractureFlowCT{Vector{Int}, Vector{Float64}}(matrix_cells, fracture_cells, connection_strength)
end

function Base.show(io::IO, d::MatrixFromFractureFlowCT)
    n = length(d.matrix_cells)
    print(io, "MatrixFromFractureFlowCT ($n connections)")
end

Jutul.symmetry(::AbstractMatrixFromFractureCT) = Jutul.CTSkewSymmetry()

function Jutul.cross_term_entities(ct::AbstractMatrixFromFractureCT, eq::ConservationLaw, model)
    return ct.matrix_cells
end

function Jutul.cross_term_entities_source(ct::AbstractMatrixFromFractureCT, eq::ConservationLaw, model)
    return ct.fracture_cells
end

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::MatrixFromFractureFlowCT, eq, dt, ldisc = local_discretization(ct, i))
    
    # Target (Matrix)
    mc = ct.matrix_cells[i]
    # Source (Fracture)
    fc = ct.fracture_cells[i]
    # Transmissibility
    T = ct.connection_strength[i]

    p_m = state_t.Pressure[mc]
    p_f = state_s.Pressure[fc]
    # Driving force: Matrix pressure - Fracture pressure
    Δp = p_m - p_f

    # Properties
    rho_m = state_t.PhaseMassDensities
    rho_f = state_s.PhaseMassDensities
    
    mob_m = state_t.PhaseMobilities
    mob_f = state_s.PhaseMobilities

    nph = size(out, 1)
    @inbounds for ph in 1:nph
        # Upwinding
        if Δp < 0
            # Fracture -> Matrix
            ρ = rho_f[ph, fc]
            λ = mob_f[ph, fc]
        else
            # Matrix -> Fracture
            ρ = rho_m[ph, mc]
            λ = mob_m[ph, mc]
        end
        out[ph] = -T * ρ * λ * Δp
    end
end

function Jutul.subcrossterm(ct::MatrixFromFractureFlowCT, ctp, m_t, m_s, map_res::Jutul.FiniteVolumeGlobalMap, map_frac::Jutul.FiniteVolumeGlobalMap, partition)
    mc = ct.matrix_cells[ctp]
    fc = ct.fracture_cells[ctp]
    T = ct.connection_strength[ctp]

    mc_local = map(c -> Jutul.local_cell(c, map_res), mc)
    fc_local = map(c -> Jutul.local_cell(c, map_frac), fc)

    return MatrixFromFractureFlowCT(mc_local, fc_local, T)
end