abstract type AbstractRelativePermeabilities <: PhaseVariables end

@enum KrScale NoKrScale TwoPointKrScale ThreePointKrScale

function Jutul.default_value(model, v::AbstractRelativePermeabilities)
    # A bit junk, but at least it's junk that sums to one for each cell.
    return 1.0/number_of_phases(model.system)
end

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

struct ReservoirRelativePermeabilities{Scaling, ph, O, OW, OG, G, R} <: AbstractRelativePermeabilities
    "Water relative permeability as a function of water saturation: ``k_{rw}(S_w)``"
    krw::O
    "Oil relative permeability (in the presence of water) as a function of oil saturation: ``k_{row}(S_o)``"
    krow::OW
    "Gas relative permeability (in the presence of water) as a function of gas saturation: ``k_{row}(S_g)``"
    krog::OG
    "Gas relative permeability as a function of gas saturation: ``k_{rg}(S_g)``"
    krg::G
    "Regions to use for each cell of the domain. Can be `nothing` if a single region is used throughout the domain."
    regions::R
    "Symbol designating the type of system, :wog for three-phase, :og for oil-gas, :wg for water-gas, etc."
    phases::Symbol
end

struct EndPointScalingCoefficients{phases} <: VectorVariables
end

function EndPointScalingCoefficients(phases::Symbol)
    return EndPointScalingCoefficients{phases}()
end

Jutul.degrees_of_freedom_per_entity(model, ::EndPointScalingCoefficients) = 4

function Jutul.default_values(model, scalers::EndPointScalingCoefficients{P}) where P
    nc = number_of_cells(model.domain)
    relperm = Jutul.get_variable(model, :RelativePermeabilities)
    kr = relperm[P]
    n = degrees_of_freedom_per_entity(model, scalers)
    kscale = zeros(n, nc)
    for i in 1:nc
        reg = JutulDarcy.region(relperm.regions, i)
        kr_i = JutulDarcy.table_by_region(kr, reg)
        (; connate, critical, s_max, k_max) = kr_i
        kscale[1, i] = connate
        kscale[2, i] = critical
        kscale[3, i] = s_max
        kscale[4, i] = k_max
    end
    return kscale
end

"""
    ReservoirRelativePermeabilities(
        w = nothing, g = nothing, ow = nothing, og = nothing,
        scaling = NoKrScale, regions = nothing)

Relative permeability with advanced features for reservoir simulation. Includes
features like rel. perm. endpoint scaling, connate water adjustment and separate
phase pair relative permeabilites for the oil phase.

# Fields

$FIELDS
"""
function ReservoirRelativePermeabilities(; w = nothing, g = nothing, ow = nothing, og = nothing, scaling = NoKrScale, regions = nothing)
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
        error("ReservoirRelativePermeabilities only implements two-phase (WO, OG, WG) or three-phase (WOG)")
    end

    F = x -> region_wrap(x, regions)
    krw = F(w)
    krow = F(ow)
    krog = F(og)
    krg = F(g)

    return ReservoirRelativePermeabilities{scaling, phases, typeof(krw), typeof(krow), typeof(krog), typeof(krg), typeof(regions)}(krw, krow, krog, krg, regions, phases)
end

function Base.getindex(m::ReservoirRelativePermeabilities, s::Symbol)
    if s == :w
        return m.krw
    elseif s == :g
        return m.krg
    elseif s == :ow
        return m.krow
    elseif s == :og
        return m.krog
    else
        error("No rel. perm. corresponding to symbol $s")
    end
end

scaling_type(::ReservoirRelativePermeabilities{T}) where T = T
scaling_type(::AbstractRelativePermeabilities) = NoKrScale

function Jutul.line_plot_data(model::SimulationModel, k::ReservoirRelativePermeabilities)
    s = collect(0:0.01:1)
    has_reg = !isnothing(k.regions)
    if has_reg
        nreg = length(k.krw)
    else
        nreg = 1
    end
    data = Matrix{Any}(undef, 2, nreg)
    ix = 1
    if !isnothing(k.krw)
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
    end
    if !isnothing(k.krg)
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
    end
    return data
end

function Jutul.subvariable(k::ReservoirRelativePermeabilities, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(k.regions, c)
    scaling = scaling_type(k)
    return ReservoirRelativePermeabilities(; w = k.krw, ow = k.krow, og = k.krog, g = k.krg, regions = regions, scaling = scaling)
end

function Base.show(io::IO, t::MIME"text/plain", kr::ReservoirRelativePermeabilities)
    println(io, "ReservoirRelativePermeabilities")
    println(io, "  functions:")
    for f in [:krw, :krow, :krog, :krg]
        k = getfield(kr, f)
        # @info k
        if isnothing(k)
            s = "(not defined)"
        else
            s = "$(length(k)) functions"
        end
        println(io, "    - $f: $s" )
    end
    if isnothing(kr.regions)
        println(io, "\n  regions: No regions defined.")
    else
        println(io, "\n  regions: $(unique(kr.regions)...).")
    end
    println(io, "\n  scaling: $(scaling_type(kr))")
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

@jutul_secondary function update_kr!(kr, relperm::ReservoirRelativePermeabilities{NoKrScale, :wog}, model, Saturations, ConnateWater, ix)
    s = Saturations
    phases = phase_indices(model.system)
    for c in ix
        @inbounds update_three_phase_relperm!(kr, relperm, phases, s, c, ConnateWater[c], nothing)
    end
    return kr
end

@jutul_secondary function update_kr!(kr, relperm::ReservoirRelativePermeabilities{NoKrScale, :wo}, model, Saturations, ix)
    s = Saturations
    regions = relperm.regions
    indices = phase_indices(model.system)
    for c in ix
        @inbounds two_phase_relperm!(kr, s, regions, relperm.krw, relperm.krow, indices, c)
    end
    return kr
end

@jutul_secondary function update_kr!(kr, relperm::ReservoirRelativePermeabilities{NoKrScale, :og}, model, Saturations, ix)
    s = Saturations
    regions = relperm.regions
    indices = phase_indices(model.system)
    for c in ix
        @inbounds two_phase_relperm!(kr, s, regions, relperm.krg, relperm.krog, reverse(indices), c)
    end
    return kr
end

@jutul_secondary function update_kr!(kr, relperm::ReservoirRelativePermeabilities{NoKrScale, :wg}, model, Saturations, ix)
    s = Saturations
    regions = relperm.regions
    indices = phase_indices(model.system)
    for c in ix
        @inbounds two_phase_relperm!(kr, s, regions, relperm.krw, relperm.krg, indices, c)
    end
    return kr
end

Base.@propagate_inbounds @inline function three_phase_oil_relperm(Krow, Krog, swcon, sg, sw)
    swc = min(swcon, value(sw) - 1e-5)
    d  = (sg + sw - swc)
    ww = (sw - swc)/d
    kro = (1-ww)*Krog + ww*Krow
    return kro
end

Base.@propagate_inbounds function two_phase_relperm!(kr, s, regions, Kr_1, Kr_2, phases, c)
    i1, i2 = phases
    reg = region(regions, c)
    sw = s[i1, c]
    sg = s[i2, c]

    kr[i1, c] = evaluate_table_by_region(Kr_1, reg, sw)
    kr[i2, c] = evaluate_table_by_region(Kr_2, reg, sg)
end

@jutul_secondary function update_kr_with_scaling!(kr, relperm::ReservoirRelativePermeabilities{<:Any, :wog}, model, Saturations, RelPermScalingW, RelPermScalingOW, RelPermScalingOG, RelPermScalingG, ConnateWater, ix)
    s = Saturations
    phases = phase_indices(model.system)
    scalers = (w = RelPermScalingW, ow = RelPermScalingOW, og = RelPermScalingOG, g = RelPermScalingG)
    for c in ix
        @inbounds update_three_phase_relperm!(kr, relperm, phases, s, c, ConnateWater[c], scalers)
    end
    return kr
end

@jutul_secondary function update_kr_with_scaling!(kr, relperm::ReservoirRelativePermeabilities{<:Any, :wo}, model, Saturations, RelPermScalingW, RelPermScalingOW, ConnateWater, ix)
    s = Saturations
    phases = phase_indices(model.system)
    scalers = (w = RelPermScalingW, ow = RelPermScalingOW)
    for c in ix
        @inbounds update_oilwater_phase_relperm!(kr, relperm, phases, s, c, ConnateWater[c], scalers)
    end
    return kr
end

Base.@propagate_inbounds @inline function update_three_phase_relperm!(kr, relperm, phase_ind, s, c, swcon, scalers)
    w, o, g = phase_ind
    reg = region(relperm.regions, c)
    krw = table_by_region(relperm.krw, reg)
    krg = table_by_region(relperm.krg, reg)
    krog = table_by_region(relperm.krog, reg)
    krow = table_by_region(relperm.krow, reg)

    sw = s[w, c]
    so = s[o, c]
    sg = s[g, c]

    Krw, Krow, Krog, Krg, swcon = three_phase_relperm(relperm, c, sw, so, sg, krw, krow, krog, krg, swcon, scalers)
    Kro = three_phase_oil_relperm(Krow, Krog, swcon, sg, sw)

    kr[w, c] = Krw
    kr[o, c] = Kro
    kr[g, c] = Krg
end

Base.@propagate_inbounds @inline function update_oilwater_phase_relperm!(kr, relperm, phase_ind, s, c, swcon, scalers)
    w, o = phase_ind
    reg = region(relperm.regions, c)
    krw = table_by_region(relperm.krw, reg)
    krow = table_by_region(relperm.krow, reg)

    sw = s[w, c]
    so = s[o, c]

    Krw, Krow = two_phase_relperm(relperm, c, sw, so, krw, krow, swcon, scalers)

    kr[w, c] = Krw
    kr[o, c] = Kro
end

function get_kr_scalers(kr::PhaseRelativePermeability)
    return (kr.connate, kr.critical, kr.s_max, kr.k_max)
end
function get_kr_scalers(scaler::AbstractMatrix, c)
    @inbounds L = scaler[1, c]
    @inbounds CR = scaler[2, c]
    @inbounds U = scaler[3, c]
    @inbounds KM = scaler[4, c]
    return (L, CR, U, KM)
end

@inline function three_phase_relperm(relperm, c, sw, so, sg, krw, krow, krog, krg, swcon, ::Nothing)
    return (krw(sw), krow(so), krog(so), krg(sg), swcon)
end

function three_phase_relperm(relperm, c, sw, so, sg, krw, krow, krog, krg, swcon, scalers)
    scaler_w, scaler_ow, scaler_og, scaler_g = scalers
    scaling = scaling_type(relperm)
    return three_phase_scaling(scaling, krw, krow, krog, krg, sw, so, sg, swcon, scaler_w, scaler_ow, scaler_og, scaler_g, c)
end

function three_phase_scaling(scaling, krw, krow, krog, krg, sw, so, sg, swcon, scaler_w, scaler_ow, scaler_og, scaler_g, c)
    L_w, CR_w, U_w, KM_w = get_kr_scalers(scaler_w, c)
    L_ow, CR_ow, U_ow, KM_ow = get_kr_scalers(scaler_ow, c)
    L_og, CR_og, U_og, KM_og = get_kr_scalers(scaler_og, c)
    L_g, CR_g, U_g, KM_g = get_kr_scalers(scaler_g, c)

    _, cr_w, u_w, km_w = get_kr_scalers(krw)
    l_w = swcon

    l_ow, cr_ow, u_ow, km_ow = get_kr_scalers(krow)
    l_og, cr_og, u_og, km_og = get_kr_scalers(krog)
    l_g, cr_g, u_g, km_g = get_kr_scalers(krg)

    R_w = 1.0 - CR_ow - L_g
    r_w = 1.0 - cr_ow - l_g

    R_g = 1.0 - CR_og - L_w
    r_g = 1.0 - cr_og - l_w

    R_ow = 1.0 - CR_w - L_g
    r_ow = 1.0 - cr_w - l_g
    U_ow = 1.0 - L_w - L_g
    u_ow = 1.0 - l_w - l_g

    R_og = 1.0 - CR_g - L_w
    r_og = 1.0 - cr_w - l_w
    U_og = 1.0 - L_g - L_w
    u_og = 1.0 - l_g - l_w

    Krw = relperm_scaling(scaling, krw, sw, cr_w, CR_w, u_w, U_w, km_w, KM_w, r_w, R_w)
    Krow = relperm_scaling(scaling, krow, so, cr_ow, CR_ow, u_ow, U_ow, km_ow, KM_ow, r_ow, R_ow)
    Krog = relperm_scaling(scaling, krog, so, cr_og, CR_og, u_og, U_og, km_og, KM_og, r_og, R_og)
    Krg = relperm_scaling(scaling, krg, sg, cr_g, CR_g, u_g, U_g, km_g, KM_g, r_g, R_g)

    return (Krw, Krow, Krog, Krg, L_w)
end

function two_phase_relperm(relperm, c, sw, so, krw, krow, swcon, scalers)
    scaler_w, scaler_ow = scalers
    scaling = scaling_type(relperm)
    return two_phase_scaling(scaling, krw, krow, sw, so, swcon, scaler_w, scaler_ow, c)
end

function two_phase_scaling(scaling, krw, krn, sw, sn, swcon, scaler_w, scaler_n, c)
    L_n, CR_n, U_n, KM_n = get_kr_scalers(scaler_n, c)
    l_n, cr_n, u_n, km_n = get_kr_scalers(krn)

    L_w, CR_w, U_w, KM_w = get_kr_scalers(scaler_w, c)
    l_w, cr_w, u_w, km_w = get_kr_scalers(krw)
    l_w = max(l_w, zero(l_w))

    R_w = 1.0 - CR_n
    r_w = 1.0 - cr_n

    R_n = 1.0 - CR_w
    r_n = 1.0 - cr_w
    U_n = 1.0 - L_w
    u_n = 1.0 - l_w

    Krw = relperm_scaling(scaling, krw, sw, cr_w, CR_w, u_w, U_w, km_w, KM_w, r_w, R_w)
    Krn = relperm_scaling(scaling, krow, sn, cr_n, CR_n, u_n, U_n, km_n, KM_n, r_n, R_n)

    return (Krw, Krn)
end

function relperm_scaling(scaling, F, s::T, cr, CR, u, U, km, KM, r, R) where T<:Real
    if scaling == ThreePointKrScale
        S = three_point_saturation_scaling(s, cr, CR, u, U, r, R)
    else
        S = two_point_saturation_scaling(s, cr, CR, u, U, r, R)
    end
    return (KM/km)*F(S)
end

function three_point_saturation_scaling(s::T, cr, CR, u, U, r, R) where T<:Real
    # @assert r >= cr
    # @assert R >= CR
    # @assert u >= r
    # @assert U >= R
    if s < CR
        S = zero(T)
    elseif s >= CR && s < R
        S = (s - CR)*(r-cr)/(R-CR) + cr
    elseif s >= R && s < U
        S = (s - R)*(u-r)/(U-R) + r
    else
        S = one(T)
    end
    return S
end

function two_point_saturation_scaling(s::T, cr, CR, u, U, r, R) where T<:Real
    # ix1 = sv < p{1};
    # ix2 = sv >= p{2};
    # ix  = ~(ix1 | ix2);
    # s_scale = (ix.*m).*s + (ix.*c + ix2);
    error("Not implemented.")
    if s < CR
        S = zero(T)
    elseif s >= CR && s < R
        S = (s - CR)*(r-cr)/(R-CR) + cr
    elseif s >= R && s < U
        S = (s - R)*(u-r)/(U-R) + r
    else
        S = one(T)
    end
    return S
end
