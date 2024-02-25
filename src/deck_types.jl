abstract type DeckPhaseVariables <: PhaseVariables end
abstract type AbstractReservoirDeckTable end
abstract type AbstractTablePVT <: AbstractReservoirDeckTable end

"""
    DeckPhaseViscosities(pvt, regions = nothing)

Secondary variable used to evaluate viscosities when a case is generated from a
input file. Typically not instantiated in user scripts.
"""
struct DeckPhaseViscosities{T, M, R} <: DeckPhaseVariables
    pvt::T
    thermal::M
    regions::R
    function DeckPhaseViscosities(pvt; regions = nothing, thermal = nothing)
        check_regions(regions)
        pvt_t = Tuple(pvt)
        if !isnothing(thermal)
            thermal = Tuple(thermal)
        end
        new{typeof(pvt_t), typeof(thermal), typeof(regions)}(pvt_t, thermal, regions)
    end
end

"""
    DeckPhaseMassDensities(pvt, regions = nothing)

Secondary variable used to evaluate densities when a case is generated from a
input file. Typically not instantiated in user scripts.
"""
struct DeckPhaseMassDensities{T, R} <: DeckPhaseVariables
    pvt::T
    regions::R
    function DeckPhaseMassDensities(pvt; regions = nothing)
        check_regions(regions)
        pvt_t = Tuple(pvt)
        new{typeof(pvt_t), typeof(regions)}(pvt_t, regions)
    end
end

"""
DeckShrinkageFactors(pvt, regions = nothing)

Secondary variable used to evaluate shrinkage factors when a case is generated
from a input file. Typically not instantiated in user scripts.
"""
struct DeckShrinkageFactors{T, W, R} <: DeckPhaseVariables
    pvt::T
    watdent::W
    regions::R
    function DeckShrinkageFactors(pvt; watdent = nothing, regions = nothing)
        check_regions(regions)
        pvt_t = Tuple(pvt)
        watdent_t = region_wrap(watdent, regions)
        new{typeof(pvt_t), typeof(watdent_t), typeof(regions)}(pvt_t, watdent, regions)
    end
end

"""
    MuBTable(pvt, regions = nothing)

Table used to evaluate viscosities and shrinkage factors when a case is
generated from a input file. Typically used to wrap tables (e.g. PVDG, PVDO) for
use in simulation.
"""
struct MuBTable{V, I}
    pressure::V
    shrinkage::V
    shrinkage_interp::I
    viscosity::V
    viscosity_interp::I
    function MuBTable(p::T, b::T, mu::T; extrapolate = true, kwarg...) where T<:AbstractVector
        @assert length(p) == length(b) == length(mu)
        I_b = get_1d_interpolator(p, b; cap_endpoints = !extrapolate, kwarg...)
        I_mu = get_1d_interpolator(p, mu; cap_endpoints = !extrapolate, kwarg...)
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

"""
    ConstMuBTable(pvtw::M) where M<:AbstractVector

Create a constant viscosity and formation-volume-factor table from a vector.
Typical usage is to wrap a PVTW type table generated from external software.
"""
function ConstMuBTable(pvtw::M) where M<:AbstractVector
    pvtw = flat_region_expand(pvtw)
    # Only one region supported atm
    pvtw = first(pvtw)
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

struct PVTO{T} <: AbstractTablePVT
    tab::T
    function PVTO(tab::T) where T<:Tuple
        for (i, t) in enumerate(tab)
            t isa PVTOTable || throw(ArgumentError("Table $i was a $(typeof(t)), should be a PVTOTable."))
        end
        return new{T}(tab)
    end
end

struct PVTOTable{T,V}
    pos::T
    rs::V
    pressure::V
    sat_pressure::V
    shrinkage::V
    viscosity::V
end

function PVTO(pvto::Vector)
    c = map(PVTOTable, pvto)
    ct = Tuple(c)
    return PVTO(ct)
end

function PVTO(pvto::PVTOTable)
    ct = (pvto, )
    return PVTO(ct)
end

function PVTO(pvto::Dict)
    return PVTO(PVTOTable(pvto))
end

function PVTOTable(d::Dict)
    rs = vec(copy(d["key"]))
    pos = vec(Int64.(d["pos"]))
    data = d["data"]
    # data, pos, rs = add_lower_pvto(data, pos, rs)
    p = vec(data[:, 1])
    B = vec(data[:, 2])
    b = 1.0 ./ B
    mu = vec(data[:, 3])
    p_sat = vec(p[pos[1:end-1]])
    T = typeof(pos)
    V = typeof(mu)
    @assert length(p) == length(b) == length(mu)
    @assert pos[end] == length(p) + 1
    @assert pos[1] == 1
    @assert length(p_sat) == length(rs) == length(pos)-1
    return PVTOTable{T, V}(pos, rs, p, p_sat, b, mu)
end

function as_printed_table(tab::PVTO, u)
    n = tab.pos[end]-1
    P = copy(tab.pressure)
    B = 1.0./tab.shrinkage
    M = copy(tab.viscosity)

    Rs = copy(tab.rs)
    # Deal with units
    InputParser.swap_unit_system!(M, u, :viscosity)
    InputParser.swap_unit_system!(Rs, u, :u_rs)
    InputParser.swap_unit_system!(P, u, :pressure)
    InputParser.swap_unit_system!(B, u, :liquid_formation_volume_factor)

    mat = Matrix{Union{Nothing, Float64}}(undef, n, 4)
    end_records = Int[]
    for i in eachindex(Rs)
        start = tab.pos[i]
        stop = (tab.pos[i+1]-1)
        pos = start:stop
        mat[start, 1] = Rs[i]
        for j in start:stop
            mat[j, 2] = P[j]
            mat[j, 3] = B[j]
            mat[j, 4] = M[j]
        end
        push!(end_records, stop)
    end
    return ("PVTO", mat, ["Rs", "Pressure", "B_o", "mu_u"], end_records)
end

function saturated_table(t::PVTOTable)
    return saturated_table(t.sat_pressure, t.rs)
end

pvt_table_vectors(pvt::PVTOTable) = (pvt.pressure, pvt.rs, pvt.sat_pressure, pvt.pos)

function shrinkage(pvt::PVTO, reg, p::T, rs, cell) where T
    tbl = table_by_region(pvt.tab, region(reg, cell))
    return interp_pvt(tbl, p, rs, tbl.shrinkage)::T
end

function viscosity(pvt::PVTO, reg, p::T, rs, cell) where T
    tbl = table_by_region(pvt.tab, region(reg, cell))
    return interp_pvt(tbl, p, rs, tbl.viscosity)::T
end

struct PVTG{T} <: AbstractTablePVT
    tab::T
    function PVTG(tab::T) where T<:Tuple
        for (i, t) in enumerate(tab)
            t isa PVTGTable || throw(ArgumentError("Table $i was a $(typeof(t)), should be a PVTGTable."))
        end
        return new{T}(tab)
    end
end

function PVTG(pvtg::Vector)
    c = map(PVTGTable, pvtg)
    ct = Tuple(c)
    return PVTG(ct)
end
# PVTG - vaporized oil
struct PVTGTable{T,V}
    pos::T
    pressure::V
    rv::V
    sat_rv::V
    shrinkage::V
    viscosity::V
end

function PVTG(pvtg::PVTGTable)
    ct = (pvtg, )
    return PVTG(ct)
end

function PVTG(pvtg::Dict)
    return PVTG(PVTGTable(pvtg))
end

function PVTGTable(d::Dict)
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
    pressure = vec(copy(d["key"]))
    # data, pos, pressure = add_lower_pvtg(data, pos, pressure)
    rv = vec(data[:, 1])
    B = vec(data[:, 2])
    b = 1.0 ./ B
    mu = vec(data[:, 3])
    rv_sat = vec(rv[pos[2:end] .- 1])
    T = typeof(pos)
    V = typeof(mu)

    @assert length(rv) == length(b) == length(mu)
    @assert pos[end] == length(rv) + 1
    @assert pos[1] == 1
    @assert length(pressure) == length(rv_sat) == length(pos)-1
    return PVTGTable{T, V}(pos, pressure, rv, rv_sat, b, mu)
end

function saturated_table(t::PVTGTable)
    return saturated_table(t.pressure, t.sat_rv)
end

function saturated_table(p, r)
    if r[1] > 0
        @assert p[1] > 0.0
        p = vcat([-1.0, 0.0], p)
        r = vcat([0.0, 0.0], r)
    end
    return get_1d_interpolator(p, r, cap_end = false)
end

pvt_table_vectors(pvt::PVTGTable) = (pvt.rv, pvt.pressure, pvt.sat_rv, pvt.pos)

function shrinkage(pvt::PVTG, reg, p::T, rv, cell) where T
    tbl = table_by_region(pvt.tab, region(reg, cell))
    # Note: Reordered arguments!
    return interp_pvt(tbl, rv, p, tbl.shrinkage)::T
end

function viscosity(pvt::PVTG, reg, p::T, rv, cell) where T
    tbl = table_by_region(pvt.tab, region(reg, cell))
    # Note: Reordered arguments!
    return interp_pvt(tbl, rv, p, tbl.viscosity)::T
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
    c = map(i -> ConstMuBTable(vec(pvtw[i, :])), axes(pvtw, 1))
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

function add_lower_pvto(data, pos, rs)
    ref_p = 101325.0
    first_offset = pos[2]-1
    start = 1:first_offset
    new_start = data[start, :]

    dp = ref_p - new_start[1, 1]
    for i in axes(new_start, 1)
        new_start[i, 1] += dp
        # new_start[i, 2] *= 0.99
    end
    @assert pos[1] == 1
    data = vcat(new_start, data)
    pos = vcat([1, first_offset+1], pos[2:end] .+ first_offset)
    rs = vcat(rs[1], rs)
    return (data, pos, rs)
end

function add_lower_pvtg(data, pos, pressure)
    ref_p = 101325.0
    first_offset = pos[2]-1
    start = 1:first_offset
    new_start = data[start, :]
    for i in axes(new_start, 1)
        # new_start[i, 2] *= 0.99
    end
    # dp = ref_p - new_start[1, 1]
    # for i in axes(new_start, 1)
        # new_start[i, 1] += dp
    # end
    @assert pos[1] == 1
    data = vcat(new_start, data)
    pos = vcat([1, first_offset+1], pos[2:end] .+ first_offset)
    pressure = vcat(ref_p, pressure)
    return (data, pos, pressure)
end

function print_deck_table!(io, tab; units = :si, self = :si, print_keyword = true)
    u_target = InputParser.DeckUnitSystem(units)
    u_self = InputParser.DeckUnitSystem(self)
    u = (from = u_self, to = u_target)
    start_label, tab_as_mat, header, end_records = as_printed_table(tab, u)

    header = copy(header)
    header[1] = "-- $(header[1])"
    function fmt(x, i, j)
        if isnothing(x)
            return ""
        else
            return Jutul.Printf.@sprintf "%4.10g" x
        end
    end
    if print_keyword && !isnothing(start_label)
        println(io, start_label)
    end
    Jutul.PrettyTables.pretty_table(
        io,
        tab_as_mat,
        tf = Jutul.PrettyTables.tf_borderless,
        hlines = :none,
        formatters = fmt,
        body_hlines_format = (' ', ' ', '/', ' '),
        body_hlines = end_records,
        alignment = :l,
        header = header
    )
    if length(end_records) > 0 && last(end_records) == size(tab_as_mat, 1)
        println(io, "    /")
    end
    println(io, "/")
end

function print_deck_table!(tab; kwarg...)
    print_deck_table!(stdout, tab; kwarg...)
end
