abstract type AbstractRelativePermeabilities <: PhaseVariables end

@enum KrScale NoKrScale TwoPointKrScale ThreePointKrScale

function Jutul.default_value(model, v::AbstractRelativePermeabilities)
    # A bit junk, but at least it's junk that sums to one for each cell.
    return 1.0/number_of_phases(model.system)
end

struct RelativePermeabilities{K, R} <: AbstractRelativePermeabilities
    relperms::K
    regions::R
    function RelativePermeabilities(kr; regions::T = nothing) where {T<:Union{Nothing, AbstractVector}}
        is_tup_tup = first(kr) isa Tuple
        check_regions(regions)
        if isnothing(regions)
            @assert !is_tup_tup || all(x -> length(x) == 1, kr) "Multiple rel. perm. functions provided, but region indicator was missing."
        end
        kr = map(x -> region_wrap(x, regions), kr)
        kr = tuple(kr...)
        return new{typeof(kr), T}(kr, regions)
    end
end

@jutul_secondary function update_kr!(kr, relperm::RelativePermeabilities, model, Saturations, ix)
    for c in ix
        reg = region(relperm.regions, c)
        @inbounds for ph in axes(kr, 1)
            s = Saturations[ph, c]
            f = table_by_region(relperm.relperms[ph], reg)
            kr[ph, c] = f(s)
        end
    end
    return kr
end

function Jutul.subvariable(k::RelativePermeabilities, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(k.regions, c)
    return RelativePermeabilities(k.relperms, regions = regions)
end

function Jutul.line_plot_data(model::SimulationModel, k::RelativePermeabilities)
    has_reg = !isnothing(k.regions)
    s = collect(0:0.01:1)
    if has_reg
        nreg = length(k.relperms[1])
    else
        nreg = 1
    end
    data = Matrix{Any}(undef, 1, nreg)
    names = phase_names(model.system)
    for r in 1:nreg
        x = []
        y = []
        labels = []
        for (i, ph) in enumerate(names)
            push!(x, s)
            push!(y, k.relperms[i][r].(s))
            push!(labels, ph)
        end
        data[1, r] = Jutul.JutulLinePlotData(x, y, labels = labels, title = "Relative permeability", xlabel = "Saturation", ylabel = "Kr")
    end
    return data
end

struct BrooksCoreyRelPerm{V, T} <: AbstractRelativePermeabilities
    exponents::V
    residuals::V
    endpoints::V
    residual_total::T
    function BrooksCoreyRelPerm(sys_or_nph::Union{MultiPhaseSystem, Integer}, exponents = 1.0, residuals = 0.0, endpoints = 1.0)
        if isa(sys_or_nph, Integer)
            nph = sys_or_nph
        else
            nph = number_of_phases(sys_or_nph)
        end
        e = expand_to_phases(exponents, nph)
        r = expand_to_phases(residuals, nph)
        epts = expand_to_phases(endpoints, nph)

        total = sum(residuals)
        new{typeof(e), typeof(total)}(e, r, epts, total)
    end
end

"""
Interpolated multiphase rel. perm. that is simple (single region, no magic for more than two phases)
"""
struct TabulatedRelPermSimple{V, M, I} <: AbstractRelativePermeabilities
    s::V
    kr::M
    interpolators::I
    function TabulatedRelPermSimple(s::AbstractVector, kr::AbstractVector; regions::Union{AbstractVector, Nothing} = nothing, kwarg...)
        nph = length(kr)
        n = length(kr[1])
        @assert nph > 0
        T = eltype(kr[1])
        #if n <= 50
        #    V = SVector{n, T}
        #else
        V = Vector{T}
        #end
        if eltype(s)<:AbstractVector
            # We got a set of different vectors that correspond to rows of kr
            @assert all(map(length, s) .== map(length, kr))
            interpolators = map((ix) -> get_1d_interpolator(V(s[ix]), V(kr[ix]); kwarg...), 1:nph)
        else
            # We got a single vector that is used for all rows
            @assert length(s) == n
            interpolators = map((ix) -> get_1d_interpolator(V(s), V(kr[ix]); kwarg...), 1:nph)
        end
        i_t = Tuple(interpolators)
        new{typeof(s), typeof(kr), typeof(i_t)}(s, kr, i_t)
    end
end

"""
Interpolated multiphase rel. perm. that is simple (single region, no magic for more than two phases)
"""
struct ReservoirRelPerm{Scaling, O, OW, OG, G, R} <: AbstractRelativePermeabilities
    krw::O
    krow::OW
    krog::OG
    krg::G
    regions::R
    phases::Symbol
end

function ReservoirRelPerm(; w = nothing, g = nothing, ow = nothing, og = nothing, scaling = NoKrScale, regions = nothing)
    has_w = !isnothing(w)
    has_g = !isnothing(g)
    has_og = !isnothing(og)
    has_ow = !isnothing(ow)
    has_o = has_og || has_ow
    if has_w && has_g && has_o
        @assert has_ow && has_og
        phases = :wog
    elseif has_w
        if has_g
            phases = :wg
        else
            @assert has_ow
            phases = :wo
        end
    elseif has_g
        if has_w
            phases = :wg
        else
            @assert has_og
            phases = :og
        end
    else
        error("ReservoirRelPerm only implements two-phase (WO, OG, WG) or three-phase (WOG)")
    end

    F = x -> region_wrap(x, regions)
    krw = F(w)
    krow = F(ow)
    krog = F(og)
    krg = F(g)

    return ReservoirRelPerm{scaling, typeof(krw), typeof(krow), typeof(krog), typeof(krg), typeof(regions)}(krw, krow, krog, krg, regions, phases)
end

function Jutul.line_plot_data(model::SimulationModel, k::ReservoirRelPerm)
    s = collect(0:0.01:1)
    has_reg = !isnothing(k.regions)
    if has_reg
        nreg = length(k.krw)
    else
        nreg = 1
    end
    data = Matrix{Any}(undef, 2, nreg)
    ix = 1
    for (krw, krow) in zip(k.krw, k.krow)
        x = []
        y = []
        labels = []
        push!(x, s)
        push!(y, krw.(s))
        push!(labels, "W")
        push!(x, 1 .- s)
        push!(y, krow.(s))
        push!(labels, "OW")
        data[1, ix] = Jutul.JutulLinePlotData(x, y, labels = labels, title = "Relative permeability", xlabel = "Water saturation", ylabel = "Water-Oil Kr")
        ix += 1
    end

    ix = 1
    for (krg, krog) in zip(k.krg, k.krog)
        x = []
        y = []
        labels = []
        push!(x, s)
        push!(y, krg.(s))
        push!(labels, "G")
        push!(x, s)
        push!(y, krog.(1 .- s))
        push!(labels, "OG")
        data[2, ix] = Jutul.JutulLinePlotData(x, y, labels = labels, title = "Relative permeability", xlabel = "Gas saturation", ylabel = "Gas-Oil Kr")
        ix += 1
    end
    return data
end

function Jutul.subvariable(k::ReservoirRelPerm, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(k.regions, c)
    swcon = Jutul.partition_variable_slice(k.swcon, c)
    return ReservoirRelPerm(k.krw, k.krow, k.krog, k.krg, swcon, regions)
end

@jutul_secondary function update_kr!(kr, kr_def::BrooksCoreyRelPerm, model, Saturations, ix)
    n, sr, kwm, sr_tot = kr_def.exponents, kr_def.residuals, kr_def.endpoints, kr_def.residual_total
    for i in ix
        for ph in axes(kr, 1)
            kr[ph, i] = brooks_corey_relperm(Saturations[ph, i], n[ph], sr[ph], kwm[ph], sr_tot)
        end
    end
    return kr
end

function brooks_corey_relperm(s::T, n::Real, sr::Real, kwm::Real, sr_tot::Real) where T
    den = 1 - sr_tot
    sat = (s - sr) / den
    sat = clamp(sat, zero(T), one(T))
    return kwm*sat^n
end

@jutul_secondary function update_kr!(kr, kr_def::TabulatedRelPermSimple, model, Saturations, ix)
    I = kr_def.interpolators
    for c in ix
        @inbounds for j in eachindex(I)
            s = Saturations[j, c]
            kr[j, c] = I[j](s)
        end
    end
    return kr
end

@jutul_secondary function update_kr!(kr, relperm::ReservoirRelPerm{NoKrScale}, model, Saturations, ix)
    s = Saturations
    phases = relperm.phases
    regions = relperm.regions
    indices = phase_indices(model.system)
    if phases == :wog
        for c in ix
            @inbounds three_phase_relperm!(kr, s, regions, relperm.krw, relperm.krg, relperm.krog, relperm.krow, indices, c)
        end
    elseif phases == :wo
        for c in ix
            @inbounds two_phase_relperm!(kr, s, regions, relperm.krw, relperm.krow, indices, c)
        end
    elseif phases == :og
        for c in ix
            @inbounds two_phase_relperm!(kr, s, regions, relperm.krg, relperm.krog, reverse(indices), c)
        end
    elseif phases == :wg
        for c in ix
            @inbounds two_phase_relperm!(kr, s, regions, relperm.krw, relperm.krg, indices, c)
        end
    end

    return kr
end

Base.@propagate_inbounds function three_phase_relperm!(kr, s, regions, Krw, Krg, Krog, Krow, phases, c)
    l, o, g = phases
    reg = region(regions, c)
    # Water
    krw = table_by_region(Krw, reg)
    sw = s[l, c]
    kr[l, c] = krw(sw)
    # Gas
    krg = table_by_region(Krg, reg)
    sg = s[g, c]
    kr[g, c] = krg(sg)
    # Oil is special
    krog = table_by_region(Krog, reg)
    krow = table_by_region(Krow, reg)
    so = s[o, c]
    swc = min(krw.connate, value(sw) - 1e-5)
    d  = (sg + sw - swc)
    ww = (sw - swc)/d
    kro = (1-ww)*krog(so) + ww*krow(so)
    kr[o, c] = kro
end

Base.@propagate_inbounds function two_phase_relperm!(kr, s, regions, Kr_1, Kr_2, phases, c)
    i1, i2 = phases
    reg = region(regions, c)
    kr1 = table_by_region(Kr_1, reg)
    sw = s[i1, c]
    kr[i1, c] = kr1(sw)
    kr2 = table_by_region(Kr_2, reg)
    sg = s[i2, c]
    kr[i2, c] = kr2(sg)
end
