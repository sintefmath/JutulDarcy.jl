struct RelativePermeabilitiesParameter <: AbstractRelativePermeabilities end

"""
    RelativePermeabilities((kr1, kr2, ...))

A simple relative permeability implementation. Assumes that each phase has a
relative permeability on the form:

``K_{r,phase} = F(S_{phase})``

Supports multiple fluid regions through the `regions` keyword.

# Examples
Single region:
```
kr1 = S -> S^2
kr2 = S -> S^3

kr = RelativePermeabilities((kr1, kr2))
```
Two regions:
```
kr1_reg1 = S -> S^2
kr2_reg1 = S -> S^3

kr1_reg2 = S -> S^3
kr2_reg2 = S -> S^4

regions # should be a vector with one entry that is 1 or 2 for each cell in the domain

kr = RelativePermeabilities(((kr1_reg1, kr2_reg1), (kr1_reg2, kr2_reg2)), regions = regions)
```
"""
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

"""
    BrooksCoreyRelativePermeabilities(
        sys_or_nph::Union{MultiPhaseSystem, Integer},
        exponents = 1.0,
        residuals = 0.0,
        endpoints = 1.0
    )

Secondary variable that implements the family of Brooks-Corey relative
permeability functions. This is a simple analytical expression for relative
permeabilities that has a limited number of parameters:

``K(S) = K_{max} * ((S - S_r)/(1 - S_r^{tot}))^N``

## Fields
$FIELDS
"""
struct BrooksCoreyRelativePermeabilities{V, T} <: AbstractRelativePermeabilities
    "Exponents for each phase"
    exponents::V
    "Residual saturations for each phase"
    residuals::V
    "Maximum relative permeability for each phase"
    endpoints::V
    "Total residual saturation over all phases"
    residual_total::T
    function BrooksCoreyRelativePermeabilities(sys_or_nph::Union{MultiPhaseSystem, Integer}, exponents = 1.0, residuals = 0.0, endpoints = 1.0)
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
    TabulatedSimpleRelativePermeabilities(s::AbstractVector, kr::AbstractVector; regions::Union{AbstractVector, Nothing} = nothing, kwarg...)

Interpolated multiphase relative permeabilities that assumes that the relative
permeability of each phase depends only on the phase saturation of that phase.
"""
struct TabulatedSimpleRelativePermeabilities{V, M, I} <: AbstractRelativePermeabilities
    s::V
    kr::M
    interpolators::I
    function TabulatedSimpleRelativePermeabilities(s::AbstractVector, kr::AbstractVector; regions::Union{AbstractVector, Nothing} = nothing, kwarg...)
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

@jutul_secondary function update_kr!(kr, kr_def::BrooksCoreyRelativePermeabilities, model, Saturations, ix)
    n, sr, kwm, sr_tot = kr_def.exponents, kr_def.residuals, kr_def.endpoints, kr_def.residual_total
    for i in ix
        for ph in axes(kr, 1)
            kr[ph, i] = brooks_corey_relperm(Saturations[ph, i], n[ph], sr[ph], kwm[ph], sr_tot)
        end
    end
    return kr
end

"""
    brooks_corey_relperm(s; n = 2.0, residual = 0.0, kr_max = 1.0, residual_total = residual)

Evaluate Brooks-Corey relative permeability function at saturation `s` for
exponent `n` and a given residual and maximum relative permeability value. If
considering a two-phase system, the total residual value over both phases should
also be passed if the other phase has a non-zero residual value.
"""
function brooks_corey_relperm(s; n = 2.0, residual = 0.0, kr_max = 1.0, residual_total = residual)
    @assert s <= 1.0
    @assert s >= 0.0
    return brooks_corey_relperm(s, n, residual, kr_max, residual_total)
end

function brooks_corey_relperm(s::T, n::Real, sr::Real, kwm::Real, sr_tot::Real) where T
    den = 1 - sr_tot
    sat = (s - sr) / den
    sat = clamp(sat, zero(T), one(T))
    return kwm*sat^n
end

@jutul_secondary function update_kr!(kr, kr_def::TabulatedSimpleRelativePermeabilities, model, Saturations, ix)
    I = kr_def.interpolators
    for c in ix
        @inbounds for j in eachindex(I)
            s = Saturations[j, c]
            kr[j, c] = I[j](s)
        end
    end
    return kr
end
