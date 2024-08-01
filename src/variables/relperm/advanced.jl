struct ReservoirRelativePermeabilities{Scaling, ph, O, OW, OG, G, R, HW, HOW, HOG, HG} <: AbstractRelativePermeabilities
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
end


"""
    ReservoirRelativePermeabilities(
        w = nothing,
        g = nothing,
        ow = nothing,
        og = nothing,
        scaling = NoKrScale(),
        regions = nothing
        hysteresis_w = NoHysteresis(),
        hysteresis_ow = NoHysteresis(),
        hysteresis_og = NoHysteresis(),
        hysteresis_g = NoHysteresis(),
        hysteresis_s_threshold = 0.0,
        hysteresis_s_eps = 1e-10
    )

Relative permeability with advanced features for reservoir simulation. Includes
features like rel. perm. endpoint scaling, connate water adjustment and separate
phase pair relative permeabilites for the oil phase.

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
        hysteresis_s_eps = 1e-10
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
        typeof(hysteresis_g)
        }(krw, krow, krog, krg, regions, phases, hysteresis_w, hysteresis_ow, hysteresis_og, hysteresis_g, scaling, hysteresis_s_threshold, hysteresis_s_eps)
end

function Jutul.get_dependencies(kr::ReservoirRelativePermeabilities, model)
    scaling = endpoint_scaling_model(kr)
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
    scaling = endpoint_scaling_model(k)
    return ReservoirRelativePermeabilities(; w = k.krw, ow = k.krow, og = k.krog, g = k.krg, regions = regions, scaling = scaling)
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
        println(io, "    - $f: $s" )
    end
    if isnothing(kr.regions)
        println(io, "\n  regions: No regions defined.")
    else
        println(io, "\n  regions: $(unique(kr.regions)...).")
    end
    println(io, "\n  scaling: $(endpoint_scaling_model(kr))")
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

        val_w = kr_hysteresis(relperm.hysteresis_w, krwd, krwi, sw, sw_max, ϵ, th)
        val_ow = kr_hysteresis(relperm.hysteresis_ow, krowd, krowi, so, so_max, ϵ, th)
        val_og = kr_hysteresis(relperm.hysteresis_og, krogd, krogi, so, so_max, ϵ, th)
        val_g = kr_hysteresis(relperm.hysteresis_g, krgd, krgi, sg, sg_max, ϵ, th)
    else
        val_w = krwd(sw)
        val_ow = krowd(so)
        val_og = krogd(so)
        val_g = krgd(sg)
    end

    kr[w, c] = val_w
    kr[o, c] = three_phase_oil_relperm(val_ow, val_og, swcon, sg, sw)
    kr[g, c] = val_g
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
    if length(f) == 1
        @assert reg == 1
        ix = 1
    else
        ix = reg + (length(f) ÷ 2)
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
    add_relperm_parameters!(model.parameters, model[:RelativePermeabilities])
    return model
end

function add_relperm_parameters!(param, kr::AbstractRelativePermeabilities)
    add_hysteresis_parameters!(param, kr)
    add_scaling_parameters!(param, kr)
    return param
end
