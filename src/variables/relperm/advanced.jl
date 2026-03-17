"""
    AbstractThreePhaseOilMethod

Abstract type for methods that combine two-phase oil relative permeabilities
``k_{row}(S_o)`` and ``k_{rog}(S_o)`` into a three-phase oil relative permeability.
"""
abstract type AbstractThreePhaseOilMethod end

"""
    SaturationWeightedOilRelperm()

Default three-phase oil relative permeability model using saturation-weighted
interpolation between ``k_{row}`` and ``k_{rog}``:

``k_{ro} = (1-w) k_{rog} + w k_{row}``

where ``w = (S_w - S_{wc})/(S_g + S_w - S_{wc})``.
"""
struct SaturationWeightedOilRelperm <: AbstractThreePhaseOilMethod end

"""
    StoneIMethod()

Stone's first model (Stone I) for three-phase oil relative permeability.
Combines two-phase oil relative permeabilities:

``k_{ro} = S_o^* \\frac{k_{row} k_{rog}}{k_{rocw}}``

where ``S_o^* = (S_o - S_{om})/(1 - S_{wc} - S_{om})`` is the normalized oil
saturation and ``k_{rocw}`` is the oil relative permeability at connate water
saturation.

Reference: Stone, H.L. "Probability Model for Estimating Three-Phase Relative
Permeability." Journal of Petroleum Technology 22.02 (1970): 214-218.
"""
struct StoneIMethod <: AbstractThreePhaseOilMethod end

"""
    StoneIIMethod()

Stone's second model (Stone II) for three-phase oil relative permeability:

``k_{ro} = k_{rocw} \\left[\\left(\\frac{k_{row}}{k_{rocw}} + k_{rw}\\right)
\\left(\\frac{k_{rog}}{k_{rocw}} + k_{rg}\\right) - (k_{rw} + k_{rg})\\right]``

Reference: Stone, H.L. "Estimation of Three-Phase Relative Permeability and
Residual Oil Data." Journal of Canadian Petroleum Technology 12.4 (1973).
"""
struct StoneIIMethod <: AbstractThreePhaseOilMethod end

struct ReservoirRelativePermeabilities{Scaling, ph, O, OW, OG, G, R, HW, HOW, HOG, HG, M} <: AbstractRelativePermeabilities
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
    "Hysteresis model for water rel. perm."
    hysteresis_w::HW
    "Hysteresis model for oil-water rel. perm."
    hysteresis_ow::HOW
    "Hysteresis model for oil-gas rel. perm."
    hysteresis_og::HOG
    "Hysteresis model for gas rel. perm."
    hysteresis_g::HG
    "Endpoint scaling model"
    scaling::Scaling
    "Threshold for hysteresis to be active"
    hysteresis_s_threshold::Float64
    "Small epsilon used in hysteresis activation check"
    hysteresis_s_eps::Float64
    "Method for computing three-phase oil relative permeability"
    three_phase_method::M
end


"""
    ReservoirRelativePermeabilities(
        w = nothing,
        g = nothing,
        ow = nothing,
        og = nothing,
        scaling = NoKrScale(),
        regions = nothing,
        hysteresis_w = NoHysteresis(),
        hysteresis_ow = NoHysteresis(),
        hysteresis_og = NoHysteresis(),
        hysteresis_g = NoHysteresis(),
        hysteresis_s_threshold = 0.0,
        hysteresis_s_eps = 1e-10,
        three_phase_method = SaturationWeightedOilRelperm()
    )

Relative permeability with advanced features for reservoir simulation. Includes
features like rel. perm. endpoint scaling, connate water adjustment and separate
phase pair relative permeabilites for the oil phase. Supports multiple methods
for combining two-phase oil relative permeabilities in three-phase systems via
`three_phase_method`:

- `SaturationWeightedOilRelperm()` (default): saturation-weighted interpolation
- `StoneIMethod()`: Stone's first model
- `StoneIIMethod()`: Stone's second model

# Fields

$FIELDS

# Examples
```julia
s = collect(range(0, 1, 100))
krw = PhaseRelativePermeability(s, s)
krog = PhaseRelativePermeability(s, s.^3)
kr_def = ReservoirRelativePermeabilities(krw = krw, krog = krog)
```
"""
function ReservoirRelativePermeabilities(;
        w = nothing,
        g = nothing,
        ow = nothing,
        og = nothing,
        scaling::AbstractKrScale = NoKrScale(),
        regions::Union{Vector{Int}, Nothing} = nothing,
        hysteresis_w::AbstractHysteresis = NoHysteresis(),
        hysteresis_ow::AbstractHysteresis = NoHysteresis(),
        hysteresis_og::AbstractHysteresis = NoHysteresis(),
        hysteresis_g::AbstractHysteresis = NoHysteresis(),
        hysteresis_s_threshold = 0.0,
        hysteresis_s_eps = 1e-10,
        three_phase_method::AbstractThreePhaseOilMethod = SaturationWeightedOilRelperm()
    )
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

    return ReservoirRelativePermeabilities{
        typeof(scaling),
        phases,
        typeof(krw),
        typeof(krow),
        typeof(krog),
        typeof(krg),
        typeof(regions),
        typeof(hysteresis_w),
        typeof(hysteresis_ow),
        typeof(hysteresis_og),
        typeof(hysteresis_g),
        typeof(three_phase_method)
        }(krw, krow, krog, krg, regions, phases, hysteresis_w, hysteresis_ow, hysteresis_og, hysteresis_g, scaling, hysteresis_s_threshold, hysteresis_s_eps, three_phase_method)
end

function Jutul.get_dependencies(kr::ReservoirRelativePermeabilities, model)
    deps = Symbol[:Saturations]
    phases = get_phases(model.system)
    has_hyst = hysteresis_is_active(kr)
    has_scaling = endpoint_scaling_is_active(kr)
    has_water = AqueousPhase() in phases
    has_oil = LiquidPhase() in phases
    has_gas = VaporPhase() in phases
    if has_hyst
        push!(deps, :MaxSaturations)
    end
    if has_water && (length(phases) > 2 || has_scaling)
        push!(deps, :ConnateWater)
    end
    if has_scaling
        if has_water
            push!(deps, :RelPermScalingW)
            if has_hyst
                push!(deps, :RelPermScalingWi)
            end
            if has_oil
                push!(deps, :RelPermScalingOW)
                if has_hyst
                    push!(deps, :RelPermScalingOWi)
                end
            end
        end
        if has_gas
            push!(deps, :RelPermScalingG)
            if has_hyst
                push!(deps, :RelPermScalingGi)
            end
            if has_oil
                push!(deps, :RelPermScalingOG)
                if has_hyst
                    push!(deps, :RelPermScalingOGi)
                end
            end
        end
    end
    out = tuple(deps...)
    return out
end

function update_secondary_variable!(kr, relperm::ReservoirRelativePermeabilities{scaling_t, ph}, model, state, ix = entity_eachindex(kr)) where {scaling_t, ph}
    s = state.Saturations
    regions = relperm.regions
    phases = phase_indices(model.system)
    scalers = get_endpoint_scalers(state, endpoint_scaling_model(relperm), Val(ph), drainage = true)
    if hysteresis_is_active(relperm)
        imb_scalers = get_endpoint_scalers(state, endpoint_scaling_model(relperm), Val(ph), drainage = false)
        s_max = state.MaxSaturations
    else
        imb_scalers = s_max = nothing
    end
    if ph == :wog
        for c in ix
            @inbounds update_three_phase_relperm!(kr, relperm, phases, s, s_max, c, state.ConnateWater[c], scalers, imb_scalers)
        end
    else
        if ph == :wo
            kr1, kr2 = relperm.krw, relperm.krow
            H1, H2 = relperm.hysteresis_w, relperm.hysteresis_ow
        elseif ph == :og
            kr1, kr2 = relperm.krog, relperm.krg
            H1, H2 = relperm.hysteresis_og, relperm.hysteresis_g
        else
            @assert ph == :wg
            kr1, kr2 = relperm.krw, relperm.krg
            H1, H2 = relperm.hysteresis_w, relperm.hysteresis_g
        end
        for c in ix
            @inbounds update_two_phase_relperm!(kr, relperm, kr1, kr2, H1, H2, phases, s, s_max, c, scalers, imb_scalers)
        end
    end
    return kr
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

function endpoint_scaling_model(x::ReservoirRelativePermeabilities)
    return x.scaling
end

function hysteresis_is_active(x::ReservoirRelativePermeabilities)
    disabled_w = x.hysteresis_w isa NoHysteresis
    disabled_ow = x.hysteresis_ow isa NoHysteresis
    disabled_og = x.hysteresis_og isa NoHysteresis
    disabled_g = x.hysteresis_g isa NoHysteresis
    return !disabled_w || !disabled_ow || !disabled_og || !disabled_g
end

function Jutul.line_plot_data(model::SimulationModel, k::ReservoirRelativePermeabilities)
    s = collect(0:0.01:1)
    has_reg = !isnothing(k.regions)
    if has_reg
        hyst = hysteresis_is_active(k)
        nreg = length(k.krw)
        if hyst
            nreg = nreg ÷ 2
        end
    else
        hyst = false
        nreg = 1
    end
    if hyst
        suffix = " (drainage)"
    else
        suffix = ""
    end
    data = Matrix{Any}(undef, 2, nreg)
    ix = 1
    if !isnothing(k.krw)
        for ix in 1:nreg
            krow = k.krow[ix]
            krw = k.krw[ix]
            x = []
            y = []
            labels = []
            push!(x, s)
            push!(y, krw.(s))
            push!(labels, "W$suffix")
            push!(x, 1 .- s)
            push!(y, krow.(s))
            push!(labels, "OW$suffix")
            if hyst
                krwi = k.krw[ix + nreg]
                krowi = k.krow[ix + nreg]
                push!(x, s)
                push!(y, krwi.(s))
                push!(labels, "W (imbibition)")
                push!(x, 1 .- s)
                push!(y, krowi.(s))
                push!(labels, "OW (imbibition)")
            end
            data[1, ix] = Jutul.JutulLinePlotData(x, y, labels = labels, title = "Relative permeability", xlabel = "Water saturation", ylabel = "Water-Oil Kr")
            ix += 1
        end
    end
    if !isnothing(k.krg)
        ix = 1
        for ix in 1:nreg
            krog = k.krog[ix]
            krg = k.krg[ix]
            x = []
            y = []
            labels = []
            push!(x, s)
            push!(y, krg.(s))
            push!(labels, "G$suffix")
            push!(x, s)
            push!(y, krog.(1 .- s))
            push!(labels, "OG$suffix")
            if hyst
                krgi = k.krg[ix + nreg]
                krogi = k.krog[ix + nreg]
                push!(x, s)
                push!(y, krgi.(s))
                push!(labels, "G (imbibition)")
                push!(x, s)
                push!(y, krogi.(1 .- s))
                push!(labels, "OG (imbibition)")
            end
            data[2, ix] = Jutul.JutulLinePlotData(x, y, labels = labels, title = "Relative permeability", xlabel = "Gas saturation", ylabel = "Gas-Oil Kr")
            ix += 1
        end
    end
    return data
end

function Jutul.subvariable(k::ReservoirRelativePermeabilities, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(k.regions, c)
    scaling = endpoint_scaling_model(k)
    return ReservoirRelativePermeabilities(; w = k.krw, ow = k.krow, og = k.krog, g = k.krg, regions = regions, scaling = scaling, three_phase_method = k.three_phase_method)
end

function Base.show(io::IO, t::MIME"text/plain", kr::ReservoirRelativePermeabilities)
    println(io, "ReservoirRelativePermeabilities")
    println(io, "  functions:")
    for f in [:w, :ow, :og, :g]
        k = getfield(kr, Symbol("kr$f"))
        if isnothing(k)
            s = "(not defined)"
        else
            hy = getfield(kr, Symbol("hysteresis_$f"))
            s = "$(length(k)) functions with $hy"
        end
        println(io, "    - kr$f: $s" )
    end
    if isnothing(kr.regions)
        println(io, "\n  regions: No regions defined.")
    else
        regstr = join(unique(kr.regions), ", ")
        println(io, "\n  regions: $regstr.")
    end
    println(io, "\n  scaling: $(endpoint_scaling_model(kr))")
    println(io, "  three-phase oil method: $(kr.three_phase_method)")
end

Base.@propagate_inbounds @inline function three_phase_oil_relperm(Krow, Krog, swcon, sg, sw, ::SaturationWeightedOilRelperm)
    swc = min(swcon, value(sw) - 1e-5)
    d  = (sg + sw - swc)
    ww = (sw - swc)/d
    kro = (1-ww)*Krog + ww*Krow
    return kro
end

Base.@propagate_inbounds @inline function three_phase_oil_relperm(Krow, Krog, swcon, sg, sw, ::StoneIMethod; krocw = 1.0, som = 0.0)
    so = 1.0 - sw - sg
    swc = swcon
    denom = 1.0 - swc - som
    if denom <= 0.0 || so <= som
        return zero(typeof(Krow))
    end
    so_star = (so - som) / denom
    if so_star <= 0
        return zero(typeof(Krow))
    end
    krocw_clamped = max(krocw, 1e-30)
    kro = so_star * Krow * Krog / krocw_clamped
    return kro
end

Base.@propagate_inbounds @inline function three_phase_oil_relperm(Krow, Krog, swcon, sg, sw, ::StoneIIMethod; krocw = 1.0, krw_val = 0.0, krg_val = 0.0)
    krocw_clamped = max(krocw, 1e-30)
    kro = krocw_clamped * ((Krow / krocw_clamped + krw_val) * (Krog / krocw_clamped + krg_val) - (krw_val + krg_val))
    kro = max(kro, zero(typeof(kro)))
    return kro
end

# Backward-compatible version that defaults to SaturationWeightedOilRelperm
Base.@propagate_inbounds @inline function three_phase_oil_relperm(Krow, Krog, swcon, sg, sw)
    return three_phase_oil_relperm(Krow, Krog, swcon, sg, sw, SaturationWeightedOilRelperm())
end

Base.@propagate_inbounds function two_phase_relperm!(kr, s, regions, Kr_1, Kr_2, phases, c)
    i1, i2 = phases
    reg = region(regions, c)
    sw = s[i1, c]
    sg = s[i2, c]

    kr[i1, c] = evaluate_table_by_region(Kr_1, reg, sw)
    kr[i2, c] = evaluate_table_by_region(Kr_2, reg, sg)
end

Base.@propagate_inbounds @inline function update_three_phase_relperm!(kr, relperm, phase_ind, s, s_max, c, swcon, scalers, scalersi)
    w, o, g = phase_ind
    reg = region(relperm.regions, c)
    krw_base = table_by_region(relperm.krw, reg)
    krg_base = table_by_region(relperm.krg, reg)
    krog_base = table_by_region(relperm.krog, reg)
    krow_base = table_by_region(relperm.krow, reg)

    if !isnothing(scalers)
        swcon = scalers.w[1, c]
    end
    krwd, krowd, krogd, krgd = get_three_phase_relperms(relperm, c, krw_base, krow_base, krog_base, krg_base, swcon, scalers)

    sw = s[w, c]
    so = s[o, c]
    sg = s[g, c]

    if hysteresis_is_active(relperm) && !isnothing(s_max)
        sw_max = s_max[w, c]
        so_max = s_max[o, c]
        sg_max = s_max[g, c]

        krwi_base = imbibition_table_by_region(relperm.krw, reg)
        krgi_base = imbibition_table_by_region(relperm.krg, reg)
        krogi_base = imbibition_table_by_region(relperm.krog, reg)
        krowi_base = imbibition_table_by_region(relperm.krow, reg)

        if isnothing(scalersi)
            swconi = swcon
        else
            swconi = scalersi.w[1, c]
        end
        krwi, krowi, krogi, krgi = get_three_phase_relperms(relperm, c, krwi_base, krowi_base, krogi_base, krgi_base, swconi, scalersi)

        ϵ = relperm.hysteresis_s_eps
        s_th = relperm.hysteresis_s_threshold
        val_w = kr_hysteresis(relperm.hysteresis_w, krwd, krwi, sw, sw_max, ϵ, s_th)
        val_ow = kr_hysteresis(relperm.hysteresis_ow, krowd, krowi, so, so_max, ϵ, s_th)
        val_og = kr_hysteresis(relperm.hysteresis_og, krogd, krogi, so, so_max, ϵ, s_th)
        val_g = kr_hysteresis(relperm.hysteresis_g, krgd, krgi, sg, sg_max, ϵ, s_th)
    else
        val_w = krwd(sw)
        val_ow = krowd(so)
        val_og = krogd(so)
        val_g = krgd(sg)
    end

    kr[w, c] = val_w
    method = relperm.three_phase_method
    kr[o, c] = compute_three_phase_oil_kr(method, val_ow, val_og, val_w, val_g, swcon, sg, sw, krowd)
    kr[g, c] = val_g
end

@inline function compute_three_phase_oil_kr(method::SaturationWeightedOilRelperm, val_ow, val_og, val_w, val_g, swcon, sg, sw, krowd)
    return three_phase_oil_relperm(val_ow, val_og, swcon, sg, sw, method)
end

@inline function compute_three_phase_oil_kr(method::StoneIMethod, val_ow, val_og, val_w, val_g, swcon, sg, sw, krowd)
    krocw = krowd.k_max
    som = krowd.critical
    return three_phase_oil_relperm(val_ow, val_og, swcon, sg, sw, method, krocw = krocw, som = som)
end

@inline function compute_three_phase_oil_kr(method::StoneIIMethod, val_ow, val_og, val_w, val_g, swcon, sg, sw, krowd)
    krocw = krowd.k_max
    return three_phase_oil_relperm(val_ow, val_og, swcon, sg, sw, method, krocw = krocw, krw_val = val_w, krg_val = val_g)
end

Base.@propagate_inbounds @inline function update_two_phase_relperm!(kr, relperm, krw, krn, H_w, H_n, phase_ind, s, s_max, c, scalers, scalersi)
    w, n = phase_ind
    reg = region(relperm.regions, c)
    krwd_base = table_by_region(krw, reg)
    krnd_base = table_by_region(krn, reg)

    krwd, krnd = get_two_phase_relperms(relperm, c, krwd_base, krnd_base, scalers)
    sw = s[w, c]
    sn = s[n, c]

    if hysteresis_is_active(relperm) && !isnothing(s_max)
        sw_max = s_max[w, c]
        sn_max = s_max[n, c]

        krwi_base = imbibition_table_by_region(krw, reg)
        krni_base = imbibition_table_by_region(krn, reg)
        krwi, krni = get_two_phase_relperms(relperm, c, krwi_base, krni_base, scalersi)

        ϵ = relperm.hysteresis_s_eps
        s_th = relperm.hysteresis_s_threshold

        val_w = kr_hysteresis(H_w, krwd, krwi, sw, sw_max, ϵ, s_th)
        val_n = kr_hysteresis(H_n, krnd, krni, sn, sn_max, ϵ, s_th)
    else
        val_w = krwd(sw)
        val_n = krnd(sn)
    end
    kr[w, c] = val_w
    kr[n, c] = val_n
end

function imbibition_table_by_region(f, reg)
    nkr = length(f)
    if nkr == 1
        @assert reg == 1
        ix = 1
    else
        # Don't go outside table range. In the case of values beyond the range,
        # some cells may use imbibiton for both curves.
        ix = min(reg + (nkr ÷ 2), nkr)
    end
    return f[ix]
end

@inline function get_three_phase_relperms(relperm, c, krw, krow, krog, krg, swcon, ::Nothing)
    return (krw, krow, krog, krg)
end

function get_three_phase_relperms(relperm, c, krw, krow, krog, krg, swcon, scalers)
    scaler_w, scaler_ow, scaler_og, scaler_g = scalers
    scaling = endpoint_scaling_model(relperm)
    return get_three_phase_scaled_relperms(scaling, krw, krow, krog, krg, swcon, scaler_w, scaler_ow, scaler_og, scaler_g, c)
end

function get_two_phase_relperms(relperm, c, krw, krow, scalers::Nothing)
    return (krw, krow)
end

function get_two_phase_relperms(relperm, c, krw, krow, scalers)
    scaler_w, scaler_ow = scalers
    scaling = endpoint_scaling_model(relperm)
    return get_two_phase_scaled_relperms(scaling, krw, krow, scaler_w, scaler_ow, c)
end

function add_relperm_parameters!(model::MultiModel)
    add_relperm_parameters!(reservoir_model(model))
    return model
end

function add_relperm_parameters!(model::SimulationModel)
    if !(model.system isa SinglePhaseSystem)
        add_relperm_parameters!(model.parameters, model[:RelativePermeabilities])
    end
    return model
end

function add_relperm_parameters!(param, kr::AbstractRelativePermeabilities)
    add_hysteresis_parameters!(param, kr)
    add_scaling_parameters!(param, kr)
    return param
end

"""
    set_relative_permeability(model;
        w = nothing,
        ow = nothing,
        og = nothing,
        g = nothing,
        swof = nothing,
        sgof = nothing,
        scaling = NoKrScale(),
        regions = nothing,
        hysteresis_w = NoHysteresis(),
        hysteresis_ow = NoHysteresis(),
        hysteresis_og = NoHysteresis(),
        hysteresis_g = NoHysteresis(),
        hysteresis_s_threshold = 0.0,
        hysteresis_s_eps = 1e-10,
        three_phase_method = SaturationWeightedOilRelperm()
    )

User-friendly helper function to set up relative permeabilities on a model. Supports
both direct specification of `PhaseRelativePermeability` objects via keyword
arguments `w`, `ow`, `og`, `g`, and also SWOF/SGOF table input.

## Table input

If `swof` is provided as a matrix with columns `[Sw, krw, krow]`, the water and
oil-water relative permeabilities are automatically constructed. Similarly, if
`sgof` is provided as a matrix with columns `[Sg, krg, krog]`, the gas and
oil-gas relative permeabilities are automatically constructed.

## Examples

Two-phase with tables:
```julia
s = collect(range(0, 1, 10))
swof = hcat(s, s.^2, reverse(s).^2)
set_relative_permeability(model, swof = swof)
```

Three-phase with Stone II:
```julia
s = collect(range(0, 1, 10))
swof = hcat(s, s.^2, reverse(s).^2)
sgof = hcat(s, s.^3, reverse(s).^3)
set_relative_permeability(model, swof = swof, sgof = sgof, three_phase_method = StoneIIMethod())
```

Direct specification:
```julia
s = collect(range(0, 1, 100))
krw = PhaseRelativePermeability(s, s)
krow = PhaseRelativePermeability(s, s.^2)
set_relative_permeability(model, w = krw, ow = krow)
```
"""
function set_relative_permeability(model;
        w = nothing,
        ow = nothing,
        og = nothing,
        g = nothing,
        swof = nothing,
        sgof = nothing,
        scaling::AbstractKrScale = NoKrScale(),
        regions::Union{Vector{Int}, Nothing} = nothing,
        hysteresis_w::AbstractHysteresis = NoHysteresis(),
        hysteresis_ow::AbstractHysteresis = NoHysteresis(),
        hysteresis_og::AbstractHysteresis = NoHysteresis(),
        hysteresis_g::AbstractHysteresis = NoHysteresis(),
        hysteresis_s_threshold = 0.0,
        hysteresis_s_eps = 1e-10,
        three_phase_method::AbstractThreePhaseOilMethod = SaturationWeightedOilRelperm()
    )
    # Convert SWOF table to relative permeability functions
    if !isnothing(swof)
        krw_t, krow_t = table_to_relperm(swof, first_label = :w, second_label = :ow)
        if isnothing(w)
            w = krw_t
        end
        if isnothing(ow)
            ow = krow_t
        end
    end

    # Convert SGOF table to relative permeability functions
    if !isnothing(sgof)
        krg_t, krog_t = table_to_relperm(sgof, first_label = :g, second_label = :og)
        if isnothing(g)
            g = krg_t
        end
        if isnothing(og)
            og = krog_t
        end
    end

    kr = ReservoirRelativePermeabilities(
        w = w,
        g = g,
        ow = ow,
        og = og,
        scaling = scaling,
        regions = regions,
        hysteresis_w = hysteresis_w,
        hysteresis_ow = hysteresis_ow,
        hysteresis_og = hysteresis_og,
        hysteresis_g = hysteresis_g,
        hysteresis_s_threshold = hysteresis_s_threshold,
        hysteresis_s_eps = hysteresis_s_eps,
        three_phase_method = three_phase_method
    )

    if model isa MultiModel
        rmodel = reservoir_model(model)
    else
        rmodel = model
    end
    set_secondary_variables!(rmodel, RelativePermeabilities = kr)
    add_relperm_parameters!(rmodel)
    return model
end
