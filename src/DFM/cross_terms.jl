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
        out[ph] = T * ρ * λ * Δp
    end
    return out
end

function Jutul.subcrossterm(ct::MatrixFromFractureFlowCT, ctp, m_t, m_s, map_res::Jutul.FiniteVolumeGlobalMap, map_frac::Jutul.FiniteVolumeGlobalMap, partition)
    mc = ct.matrix_cells[ctp]
    fc = ct.fracture_cells[ctp]
    T = ct.connection_strength[ctp]

    mc_local = map(c -> Jutul.local_cell(c, map_res), mc)
    fc_local = map(c -> Jutul.local_cell(c, map_frac), fc)

    return MatrixFromFractureFlowCT(mc_local, fc_local, T)
end

struct FracturesFromWellFlowCT{I<:AbstractVector} <: AbstractReservoirFromWellCT
    fracture_cells::I
    well_cells::I
end

function Jutul.cross_term_entities(ct::FracturesFromWellFlowCT, eq::ConservationLaw, model)
    return ct.fracture_cells
end

# Jutul.can_impact_cross_term(force_t::PerforationMask, cross_term::AbstractReservoirFromWellCT) = true

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::FracturesFromWellFlowCT, eq, dt, ldisc = local_discretization(ct, i))
    sys = model_t.system
    rhoS = reference_densities(sys)
    conn = cross_term_perforation_get_fracture_conn(ct, i, state_s, state_t)
    # Call smaller interface that is easy to specialize
    if haskey(state_s, :MassFractions)
        @inbounds simple_well_perforation_flux!(out, sys, state_t, state_s, rhoS, conn)
    else
        @inbounds multisegment_well_perforation_flux!(out, sys, state_t, state_s, rhoS, conn)
    end
    return out
end

function cross_term_perforation_get_fracture_conn(ct, i, state_s, state_t)
    @inbounds begin 
        fracture_cell = ct.fracture_cells[i]
        well_cell = ct.well_cells[i]
        WI = state_s.FractureWellIndices[i]
        gdz = state_s.PerforationGravityDifference[i]
        p_well = state_s.Pressure
        p_res = state_t.Pressure
        dp = p_well[well_cell] - p_res[fracture_cell]
    end

    # Wrap the key connection data in tuple for easy extension later
    conn = (
        dp = dp,
        WI = WI,
        gdz = gdz,
        well = well_cell,
        perforation = i,
        reservoir = fracture_cell
    )
    return conn
end