struct ReservoirRelativePermeabilities{Scaling, ph, O, OW, OG, G, R, H} <: AbstractRelativePermeabilities
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
    "Hysteresis model for each phase"
    hysteresis::H
    "Endpoint scaling model"
    scaling::Scaling
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
        hysteresis_o = NoHysteresis(),
        hysteresis_g = NoHysteresis()
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
        hysteresis_o::AbstractHysteresis = NoHysteresis(),
        hysteresis_g::AbstractHysteresis = NoHysteresis()
    )
    has_w = !isnothing(w)
    has_g = !isnothing(g)
    has_og = !isnothing(og)
    has_ow = !isnothing(ow)
    has_o = has_og || has_ow
    if has_w && has_g && has_o
        @assert has_ow && has_og
        phases = :wog
        hyst = (hysteresis_w, hysteresis_o, hysteresis_g)
    elseif has_w
        if has_g
            phases = :wg
            hyst = (hysteresis_w, hysteresis_g)
        else
            @assert has_ow
            phases = :wo
            hyst = (hysteresis_w, hysteresis_o)
        end
    elseif has_g
        if has_w
            phases = :wg
            hyst = (hysteresis_w, hysteresis_g)
        else
            @assert has_og
            phases = :og
            hyst = (hysteresis_o, hysteresis_g)
        end
    else
        error("ReservoirRelativePermeabilities only implements two-phase (WO, OG, WG) or three-phase (WOG)")
    end

    F = x -> region_wrap(x, regions)
    krw = F(w)
    krow = F(ow)
    krog = F(og)
    krg = F(g)

    return ReservoirRelativePermeabilities{typeof(scaling), phases, typeof(krw), typeof(krow), typeof(krog), typeof(krg), typeof(regions), typeof(hyst)}(krw, krow, krog, krg, regions, phases, hyst, scaling)
end

function Jutul.get_dependencies(kr::ReservoirRelativePermeabilities, model)
    scaling = endpoint_scaling_model(kr)
    deps = Symbol[:Saturations]
    phases = get_phases(model.system)
    has_scaling = endpoint_scaling_is_active(kr)
    has_water = AqueousPhase() in phases
    has_oil = LiquidPhase() in phases
    has_gas = VaporPhase() in phases

    if has_water && (length(phases) > 2 || has_scaling)
        push!(deps, :ConnateWater)
    end
    if has_scaling
        if has_water
            push!(deps, :RelPermScalingW)
            if has_oil
                push!(deps, :RelPermScalingOW)
            end
        end
        if has_gas
            push!(deps, :RelPermScalingG)
            if has_oil
                push!(deps, :RelPermScalingOG)
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
    if !endpoint_scaling_is_active(relperm)
        if ph == :wog
            for c in ix
                @inbounds update_three_phase_relperm!(kr, relperm, phases, s, c, state.ConnateWater[c], nothing)
            end
        elseif ph == :wo
            for c in ix
                @inbounds two_phase_relperm!(kr, s, regions, relperm.krw, relperm.krow, phases, c)
            end
        elseif ph == :og
            for c in ix
                @inbounds two_phase_relperm!(kr, s, regions, relperm.krg, relperm.krog, reverse(phases), c)
            end
        else
            @assert ph == :wg
            for c in ix
                @inbounds two_phase_relperm!(kr, s, regions, relperm.krw, relperm.krg, phases, c)
            end
        end
    else
        if ph == :wog
            scalers = (
                w = state.RelPermScalingW,
                ow = state.RelPermScalingOW,
                og = state.RelPermScalingOG,
                g = state.RelPermScalingG
            )
            for c in ix
                @inbounds update_three_phase_relperm!(kr, relperm, phases, s, c, state.ConnateWater[c], scalers)
            end
        elseif ph == :wo
            scalers = (
                w = state.RelPermScalingW,
                ow = state.RelPermScalingOW
            )
            for c in ix
                @inbounds update_oilwater_phase_relperm!(kr, relperm, phases, s, c, state.ConnateWater[c], scalers)
            end
        elseif ph == :og
            scalers = (
                og = state.RelPermScalingOG,
                g = state.RelPermScalingG
                )
            for c in ix
                @inbounds update_oilgas_phase_relperm!(kr, relperm, phases, s, c, scalers)
            end
        else
            @assert ph == :wg
            scalers = (
                w = state.RelPermScalingW,
                g = state.RelPermScalingG
            )
            for c in ix
                @inbounds update_gaswater_phase_relperm!(kr, relperm, phases, s, c, state.ConnateWater[c], scalers)
            end
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
    for f in [:krw, :krow, :krog, :krg]
        k = getfield(kr, f)
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
    kr[o, c] = Krow
end

Base.@propagate_inbounds @inline function update_gaswater_phase_relperm!(kr, relperm, phase_ind, s, c, swcon, scalers)
    w, g = phase_ind
    reg = region(relperm.regions, c)
    krw = table_by_region(relperm.krw, reg)
    krg = table_by_region(relperm.krg, reg)

    sw = s[w, c]
    sg = s[g, c]

    Krw, Krg = two_phase_relperm(relperm, c, sw, sg, krw, krg, swcon, scalers)

    kr[w, c] = Krw
    kr[g, c] = Krg
end

Base.@propagate_inbounds @inline function update_oilgas_phase_relperm!(kr, relperm, phase_ind, s, c, scalers)
    o, g = phase_ind
    reg = region(relperm.regions, c)
    krog = table_by_region(relperm.krog, reg)
    krg = table_by_region(relperm.krg, reg)

    so = s[o, c]
    sg = s[g, c]

    Krog, Krg = two_phase_relperm(relperm, c, so, sg, krog, krg, missing, scalers)

    kr[o, c] = Krog
    kr[g, c] = Krg
end


@inline function three_phase_relperm(relperm, c, sw, so, sg, krw, krow, krog, krg, swcon, ::Nothing)
    return (krw(sw), krow(so), krog(so), krg(sg), swcon)
end

function three_phase_relperm(relperm, c, sw, so, sg, krw, krow, krog, krg, swcon, scalers)
    scaler_w, scaler_ow, scaler_og, scaler_g = scalers
    scaling = endpoint_scaling_model(relperm)
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

    # Residual water
    R_w = 1.0 - CR_ow - L_g
    r_w = 1.0 - cr_ow - l_g
    # Residual gas
    R_g = 1.0 - CR_og - L_w
    r_g = 1.0 - cr_og - l_w
    # Residual oil
    R_ow = 1.0 - CR_w - L_g
    r_ow = 1.0 - cr_w - l_g
    R_og = 1.0 - CR_g - L_w
    r_og = 1.0 - cr_w - l_w
    # Maximum saturation oil (in persence of water)
    U_ow = 1.0 - L_w - L_g
    u_ow = 1.0 - l_w - l_g
    # Maximum saturation oil (in persence of gas)
    U_og = 1.0 - L_g - L_w
    u_og = 1.0 - l_g - l_w

    Krw = relperm_scaling(scaling, krw, sw, cr_w, CR_w, u_w, U_w, km_w, KM_w, r_w, R_w)
    Krow = relperm_scaling(scaling, krow, so, cr_ow, CR_ow, u_ow, U_ow, km_ow, KM_ow, r_ow, R_ow)
    Krog = relperm_scaling(scaling, krog, so, cr_og, CR_og, u_og, U_og, km_og, KM_og, r_og, R_og)
    Krg = relperm_scaling(scaling, krg, sg, cr_g, CR_g, u_g, U_g, km_g, KM_g, r_g, R_g)


    ScaledPhaseRelativePermeability(krw, scaling, connate = L_w, critical = CR_w, k_max = KM_w, s_max = U_w)
    return (Krw, Krow, Krog, Krg, L_w)
end

function two_phase_relperm(relperm, c, sw, so, krw, krow, swcon, scalers)
    scaler_w, scaler_ow = scalers
    scaling = endpoint_scaling_model(relperm)
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
    Krn = relperm_scaling(scaling, krn, sn, cr_n, CR_n, u_n, U_n, km_n, KM_n, r_n, R_n)

    return (Krw, Krn)
end

