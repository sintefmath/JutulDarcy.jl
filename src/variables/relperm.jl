abstract type AbstractRelativePermeabilities <: PhaseVariables end

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
struct ThreePhaseRelPerm{O, OW, OG, G, S, R} <: AbstractRelativePermeabilities
    krw::O
    krow::OW
    krog::OG
    krg::G
    swcon::S
    regions::R
end

function ThreePhaseRelPerm(; w, g, ow, og, swcon = 0.0, regions = nothing)
    F = x -> region_wrap(x, regions)
    return ThreePhaseRelPerm(F(w), F(ow), F(og), F(g), swcon, regions)
end

function Jutul.line_plot_data(model::SimulationModel, k::ThreePhaseRelPerm)
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

function Jutul.subvariable(k::ThreePhaseRelPerm, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(k.regions, c)
    swcon = Jutul.partition_variable_slice(k.swcon, c)
    return ThreePhaseRelPerm(k.krw, k.krow, k.krog, k.krg, swcon, regions)
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

@jutul_secondary function update_kr!(kr, relperm::ThreePhaseRelPerm, model, Saturations, ix)
    s = Saturations
    swcon = relperm.swcon

    l, o, g = phase_indices(model.system)
    @inbounds for c in ix
        reg = region(relperm.regions, c)
        # Water
        krw = table_by_region(relperm.krw, reg)
        sw = s[l, c]
        kr[l, c] = krw(sw)
        # Gas
        krg = table_by_region(relperm.krg, reg)
        sg = s[g, c]
        kr[g, c] = krg(sg)
        # Oil is special
        krog = table_by_region(relperm.krog, reg)
        krow = table_by_region(relperm.krow, reg)
        so = s[o, c]
        swc = min(swcon[reg], value(sw) - 1e-5)
        d  = (sg + sw - swc)
        ww = (sw - swc)/d
        kro = (1-ww)*krog(so) + ww*krow(so)
        kr[o, c] = kro
    end
    return kr
end
