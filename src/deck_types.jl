abstract type DeckPhaseVariables <: PhaseVariables end
abstract type AbstractReservoirDeckTable end

export MuBTable, ConstMuBTable

abstract type AbstractTablePVT <: AbstractReservoirDeckTable end


struct DeckViscosity{T, R} <: DeckPhaseVariables
    pvt::T
    regions::R
    function DeckViscosity(pvt; regions = nothing)
        check_regions(regions)
        pvt_t = Tuple(pvt)
        new{typeof(pvt_t), typeof(regions)}(pvt_t, regions)
    end
end

export DeckDensity, RelativePermeabilities, ThreePhaseCompositionalDensitiesLV, PhaseMassFractions, PhaseMassFractions
export ThreePhaseLBCViscositiesLV
struct DeckDensity{T, R} <: DeckPhaseVariables
    pvt::T
    regions::R
    function DeckDensity(pvt; regions = nothing)
        check_regions(regions)
        pvt_t = Tuple(pvt)
        new{typeof(pvt_t), typeof(regions)}(pvt_t, regions)
    end
end

struct DeckShrinkageFactors{T, R} <: DeckPhaseVariables
    pvt::T
    regions::R
    function DeckShrinkageFactors(pvt; regions = nothing)
        check_regions(regions)
        pvt_t = Tuple(pvt)
        new{typeof(pvt_t), typeof(regions)}(pvt_t, regions)
    end
end


struct MuBTable{V, I}
    pressure::V
    shrinkage::V
    shrinkage_interp::I
    viscosity::V
    viscosity_interp::I
    function MuBTable(p::T, b::T, mu::T; kwarg...) where T<:AbstractVector
        @assert length(p) == length(b) == length(mu)
        I_b = get_1d_interpolator(p, b; kwarg...)
        I_mu = get_1d_interpolator(p, mu; kwarg...)
        new{T, typeof(I_b)}(p, b, I_b, mu, I_mu)
    end
end

function MuBTable(pvtx::T; kwarg...) where T<:AbstractMatrix
    N = size(pvtx, 1)
    p = vec(pvtx[:, 1])
    B = vec(pvtx[:, 2])
    b = 1.0./B
    mu = vec(pvtx[:, 3])

    # V = SVector{N, eltype(mu)}
    V = Vector{eltype(mu)}
    MuBTable(V(p), V(b), V(mu); kwarg...)
end

function viscosity(tbl::MuBTable, p)
    return tbl.viscosity_interp(p)
end

function shrinkage(tbl::MuBTable, p)
    return tbl.shrinkage_interp(p)
end

struct ConstMuBTable{R}
    p_ref::R
    b_ref::R
    b_c::R
    mu_ref::R
    mu_c::R
end

function ConstMuBTable(pvtw::M) where M<:AbstractVector
    return ConstMuBTable(pvtw[1], 1.0/pvtw[2], pvtw[3], pvtw[4], pvtw[5])
end

function viscosity(pvt::AbstractTablePVT, reg, p, cell)
    tbl = table_by_region(pvt.tab, region(reg, cell))
    return viscosity(tbl, p)
end


function viscosity(tbl::ConstMuBTable, p::T) where T
    p_r = tbl.p_ref
    μ_r = tbl.mu_ref
    c = tbl.mu_c

    F = -c*(p - p_r)
    μ = μ_r/(one(T) + F + 0.5*F^2)
    return μ::T
end

# 
function shrinkage(pvt::AbstractTablePVT, reg, p::T, cell) where T
    tbl = table_by_region(pvt.tab, region(reg, cell))
    return shrinkage(tbl, p)::T
end

function shrinkage(tbl::ConstMuBTable, p::T) where T
    p_r = tbl.p_ref
    b_r = tbl.b_ref
    c = tbl.b_c

    F = c*(p - p_r)
    b = b_r*(one(T) + F + 0.5*F^2)
    return b::T
end

# PVTO - dissolved gas
struct PVTO{T,V} <: AbstractTablePVT where {T<:AbstractArray, V<:AbstractArray}
    pos::T
    rs::V
    pressure::V
    sat_pressure::V
    shrinkage::V
    viscosity::V
end

function PVTO(d::Dict)
    rs = vec(copy(d["key"]))
    pos = vec(Int64.(d["pos"]))
    data = d["data"]
    p = vec(data[:, 1])
    B = vec(data[:, 2])
    b = 1.0 ./ B
    mu = vec(data[:, 3])
    p_sat = vec(p[pos[1:end-1]])
    T = typeof(pos)
    V = typeof(mu)
    return PVTO{T, V}(pos, rs, p, p_sat, b, mu)
end

pvt_table_vectors(pvt::PVTO) = (pvt.pressure, pvt.rs, pvt.sat_pressure, pvt.pos)

function shrinkage(pvt::PVTO, reg, p::T, rs, cell) where T
    return interp_pvt(pvt, p, rs, :shrinkage)::T
end

function viscosity(pvt::PVTO, reg, p::T, rs, cell) where T
    return interp_pvt(pvt, p, rs, :viscosity)::T
end


# PVTG - vaporized oil
struct PVTG{T,V} <: AbstractTablePVT where {T<:AbstractArray, V<:AbstractArray}
    pos::T
    pressure::V
    rv::V
    sat_rv::V
    shrinkage::V
    viscosity::V
end

function PVTG(d::Dict)
    pressure = vec(copy(d["key"]))
    pos = vec(Int64.(d["pos"]))
    data = copy(d["data"])
    for i in 1:length(pos)-1
        start = pos[i]
        stop = pos[i+1]-1
        if stop - start > 0
            if data[start, 1] > data[start+1, 1]
                # Reverse table 
                data[start:stop, :] = data[stop:-1:start, :]
            end
        end
    end
    rv = vec(data[:, 1])
    B = vec(data[:, 2])
    b = 1.0 ./ B
    mu = vec(data[:, 3])
    rv_sat = vec(rv[pos[2:end] .- 1])
    T = typeof(pos)
    V = typeof(mu)
    return PVTG{T, V}(pos, pressure, rv, rv_sat, b, mu)
end

pvt_table_vectors(pvt::PVTG) = (pvt.rv, pvt.pressure, pvt.sat_rv, pvt.pos)

function shrinkage(pvt::PVTG, reg, p::T, rv, cell) where T
    # Note: Reordered arguments!
    return interp_pvt(pvt, rv, p, :shrinkage)::T
end

function viscosity(pvt::PVTG, reg, p::T, rv, cell) where T
    # Note: Reordered arguments!
    return interp_pvt(pvt, rv, p, :viscosity)::T
end

struct PVDO{T} <: AbstractTablePVT
    tab::T
end

function PVDO(pvdo::AbstractArray)
    c = map(MuBTable, pvdo)
    ct = Tuple(c)
    PVDO{typeof(ct)}(ct)
end

struct PVDG{T} <: AbstractTablePVT
    tab::T
end

function PVDG(pvdo::AbstractArray)
    c = map(MuBTable, pvdo)
    ct = Tuple(c)
    PVDG{typeof(ct)}(ct)
end

struct PVTW_EXTENDED{T} <: AbstractTablePVT
    tab::T
end

function PVTW_EXTENDED(pvtw_extended::AbstractArray)
    c = map(MuBTable, pvtw_extended)
    ct = Tuple(c)
    PVTW_EXTENDED{typeof(ct)}(ct)
end

struct PVTW{N, T} <: AbstractTablePVT
    tab::NTuple{N, T}
end

function PVTW(pvtw::AbstractArray)
    if eltype(pvtw)<:AbstractFloat
        pvtw = [pvtw]
    end
    c = map(x -> ConstMuBTable(vec(x)), pvtw)
    ct = Tuple(c)
    N = length(ct)
    T = typeof(ct[1])
    PVTW{N, T}(ct)
end

struct PVCDO{N, T} <: AbstractTablePVT
    tab::NTuple{N, T}
end

function PVCDO(pvcdo::AbstractArray)
    if eltype(pvcdo)<:AbstractFloat
        pvcdo = [pvcdo]
    end
    c = map(x -> ConstMuBTable(vec(x)), pvcdo)
    ct = Tuple(c)
    N = length(c)
    N = length(ct)
    T = typeof(ct[1])
    PVCDO{N, T}(ct)
end


struct LinearlyCompressiblePoreVolume{R} <: ScalarVariable where {R<:Real}
    reference_pressure::R
    expansion::R
    function LinearlyCompressiblePoreVolume(; reference_pressure::T = 101325.0, expansion::T = 1e-10) where T
        new{T}(reference_pressure, expansion)
    end
end


# abstract type AbstractTableSaturation <: AbstractTableDeck end

# struct RelativePermeabilityTable <: AbstractTableSaturation

# end

# struct CapillaryPressureTable <: AbstractTableSaturation

# end

# Regions go in the outer part

