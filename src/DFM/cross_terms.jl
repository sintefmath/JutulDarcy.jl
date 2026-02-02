abstract type AbstractMatrixFromFractureCT <: Jutul.AdditiveCrossTerm end

Jutul.symmetry(::AbstractMatrixFromFractureCT) = Jutul.CTSkewSymmetry()

function Jutul.cross_term_entities(ct::AbstractMatrixFromFractureCT, eq::ConservationLaw, model)
    return ct.matrix_cells
end

function Jutul.cross_term_entities_source(ct::AbstractMatrixFromFractureCT, eq::ConservationLaw, model)
    return ct.fracture_cells
end
struct MatrixFromFractureFlowCT{I<:AbstractVector, T<:AbstractVector} <: AbstractMatrixFromFractureCT
    matrix_cells::I
    fracture_cells::I
    connection_strength_flow::T
    gravity_potential::T
end

function MatrixFromFractureFlowCT(matrix_cells::Vector{Int}, fracture_cells::Vector{Int}, connection_strength_flow::Vector{Float64}, gravity_potential::Vector{Float64})
    return MatrixFromFractureFlowCT{Vector{Int}, Vector{Float64}}(matrix_cells, fracture_cells, connection_strength_flow, gravity_potential)
end

function Base.show(io::IO, d::MatrixFromFractureFlowCT)
    n = length(d.matrix_cells)
    print(io, "MatrixFromFractureFlowCT ($n connections)")
end

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::MatrixFromFractureFlowCT, eq, dt, ldisc = local_discretization(ct, i))
    
    conn = cross_term_matrix_fracture_get_conn(ct, i, state_s, state_t)
    sys = model_t.system
    nph = number_of_phases(sys)
    @inbounds begin 
        for ph in 1:nph
            out[ph] = matrix_fracture_phase_mass_flux(conn, state_t, state_s, ph)
        end
    end
    return out
end

function cross_term_matrix_fracture_get_conn(ct, i, state_f, state_m)
    @inbounds begin 
        matrix_cell = ct.matrix_cells[i]
        fracture_cell = ct.fracture_cells[i]
        T = ct.connection_strength_flow[i]
        if hasproperty(ct, :connection_strength_thermal)
            T_th = ct.connection_strength_thermal[i]
        else
            T_th = missing
        end
        gdz = ct.gravity_potential[i]
        p_m = state_m.Pressure
        p_f = state_f.Pressure
        dp = p_f[fracture_cell] - p_m[matrix_cell]
    end
    # Wrap the key connection data in tuple for easy extension later
    conn = (
        dp = dp,
        T = T,
        T_th = T_th,
        gdz = gdz,
        matrix = matrix_cell,
        fracture = fracture_cell
    )
    return conn
end

function matrix_fracture_phase_mass_flux(conn, state_m, state_f, ph)
    # ψ is pressure difference from reservoir to well. If it is negative, we are injecting into the reservoir.
    ψ = matrix_fracture_phase_potential_difference(conn, state_m, state_f, ph)
    if ψ < 0
        cf = conn.fracture
        # Fracture -> Matrix
        if haskey(state_f, :PhaseMassMobilities)
            ρλ = state_f.PhaseMassMobilities[ph, cf]
        else
            ρ = state_f.PhaseMassDensities[ph, cf]
            λ = state_f.PhaseMobilities[ph, cf]
            ρλ = ρ*λ
        end
    else
        cm = conn.matrix
        # Matrix -> Fracture
        if haskey(state_m, :PhaseMassMobilities)
            ρλ = state_m.PhaseMassMobilities[ph, cm]
        else
            ρ = state_m.PhaseMassDensities[ph, cm]
            λ = state_m.PhaseMobilities[ph, cm]
            ρλ = ρ*λ
        end
    end
    return ρλ*ψ
end

function matrix_fracture_phase_potential_difference(conn, state_m, state_f, ix)
    dp = conn.dp
    T = conn.T
    T, dp = Base.promote(T, dp)
    if haskey(state_m, :PermeabilityMultiplier)
        K_mul = state_m[:PermeabilityMultiplier][conn.matrix]
        T *= K_mul
    end
    if conn.gdz != 0.0
        if haskey(state_f, :ConnectionPressureDrop)
            dp += state_f.ConnectionPressureDrop[conn.fracture]
        else
            ρ_m = state_m.PhaseMassDensities[ix, conn.matrix]
            if haskey(state_f, :PhaseMassDensities)
                ρ_f = state_f.PhaseMassDensities[ix, conn.fracture]
                ρ = 0.5*(ρ_m + ρ_f)
            else
                ρ = ρ_m
            end
            dp += ρ*conn.gdz
        end
    end
    return -T*dp
end

function Jutul.subcrossterm(ct::MatrixFromFractureFlowCT, ctp, m_t, m_s, map_res::Jutul.FiniteVolumeGlobalMap, map_frac::Jutul.FiniteVolumeGlobalMap, partition)
    mc = ct.matrix_cells[ctp]
    fc = ct.fracture_cells[ctp]
    T = ct.connection_strength[ctp]

    mc_local = map(c -> Jutul.local_cell(c, map_res), mc)
    fc_local = map(c -> Jutul.local_cell(c, map_frac), fc)

    return MatrixFromFractureFlowCT(mc_local, fc_local, T)
end

struct MatrixFromFractureThermalCT{I<:AbstractVector, T<:AbstractVector} <: AbstractMatrixFromFractureCT
    matrix_cells::I
    fracture_cells::I
    connection_strength_flow::T
    connection_strength_thermal::T
    gravity_potential::T
end

function MatrixFromFractureThermalCT(matrix_cells::Vector{Int}, fracture_cells::Vector{Int}, connection_strength_flow::Vector{Float64}, connection_strength_thermal::Vector{Float64}, gravity_potential::Vector{Float64})
    return MatrixFromFractureThermalCT{Vector{Int}, Vector{Float64}}(matrix_cells, fracture_cells, connection_strength_flow, connection_strength_thermal, gravity_potential)
end

function Base.show(io::IO, d::MatrixFromFractureThermalCT)
    n = length(d.matrix_cells)
    print(io, "MatrixFromFractureThermalCT ($n connections)")
end


function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::MatrixFromFractureThermalCT, eq, dt, ldisc = local_discretization(ct, i))
    
    # Target (Matrix)
    mc = ct.matrix_cells[i]
    # Source (Fracture)
    fc = ct.fracture_cells[i]
    # Transmissibility
    Λ = ct.connection_strength_thermal[i]

    T_m = state_t.Temperature[mc]
    T_f = state_s.Temperature[fc]
    # Driving force: Matrix pressure - Fracture pressure
    ΔT = T_m - T_f

    # Properties
    h_m = state_t.FluidEnthalpy
    h_f = state_s.FluidEnthalpy
    
    nph = size(out, 1)
    @inbounds for ph in 1:nph
        # Upwinding
        if ΔT < 0
            # Fracture -> Matrix
            h = h_f[ph, fc]
        else
            # Matrix -> Fracture
            h = h_m[ph, mc]
        end

        out[ph] = Λ*ΔT
    end
    return out
end

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::MatrixFromFractureThermalCT, eq, dt, ldisc = local_discretization(ct, i))
    
    conn = cross_term_matrix_fracture_get_conn(ct, i, state_s, state_t)
    sys = model_t.system
    nph = number_of_phases(sys)
    qh = matrix_fracture_thermal_flux(conn, state_t, state_s, nph)
    out[] = qh
    return out
end


function matrix_fracture_thermal_flux(conn, state_m, state_f, nph)
    
    # Get connection properties
    cf = conn.fracture
    cm = conn.matrix
    Λ = conn.T_th
    # Compute conductive heat flux
    T_m = state_m.Temperature[cm]
    T_f = state_f.Temperature[cf]
    ΔT = T_f - T_m
    qh = -Λ*ΔT
    # Add advective heat flux
    h_m = state_m.FluidEnthalpy
    h_f = state_f.FluidEnthalpy
    @inbounds for ph in 1:nph
        qm = matrix_fracture_phase_mass_flux(conn, state_m, state_f, ph)
        if qm < 0
            # Fracture -> Matrix
            h = h_f[ph, cf]
        else
            # Matrix -> Fracture
            h = h_m[ph, cm]
        end
        qh += qm*h
    end

    return qh
    
end

function Jutul.subcrossterm(ct::MatrixFromFractureThermalCT, ctp, m_t, m_s, map_res::Jutul.FiniteVolumeGlobalMap, map_frac::Jutul.FiniteVolumeGlobalMap, partition)
    mc = ct.matrix_cells[ctp]
    fc = ct.fracture_cells[ctp]
    T = ct.connection_strength[ctp]

    mc_local = map(c -> Jutul.local_cell(c, map_res), mc)
    fc_local = map(c -> Jutul.local_cell(c, map_frac), fc)

    return MatrixFromFractureThermalCT(mc_local, fc_local, T)
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
        error("FractureFromWellFlowCT currently only supports multisegment wells.")
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

struct FracturesFromWellThermalCT{I<:AbstractVector} <: AbstractReservoirFromWellCT
    fracture_cells::I
    well_cells::I
end

function Jutul.cross_term_entities(ct::FracturesFromWellThermalCT, eq::ConservationLaw, model)
    return ct.fracture_cells
end

function update_cross_term_in_entity!(out, i,
    state_res, state0_res,
    state_well, state0_well, 
    model_res, model_well,
    ct::FracturesFromWellThermalCT, eq, dt, ldisc = local_discretization(ct, i))
    # Unpack properties
    sys = model_res.system
    nph = number_of_phases(sys)
    @inbounds begin 
        fracture_cell = ct.fracture_cells[i]
        well_cell = ct.well_cells[i]
        WI = state_well.FractureWellIndices[i]
        WIth = state_well.FractureWellThermalIndices[i]
        gdz = state_well.PerforationGravityDifference[i]
        p_well = state_well.Pressure
        p_res = state_res.Pressure
        dp = p_well[well_cell] - p_res[fracture_cell]
        conn = (
            dp = dp,
            WI = WI,
            WIth = WIth,
            gdz = gdz,
            well = well_cell,
            perforation = i,
            reservoir = fracture_cell
        )
    end

    λ_t = sum(perforation_reservoir_mobilities(state_res, state_well, sys, reservoir_cell, well_cell))
    qh = perforation_phase_thermal_flux(λ_t, conn, state_res, state_well, nph)
    out[] = qh

end
