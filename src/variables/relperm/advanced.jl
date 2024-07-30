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
    scalers = get_endpoint_scalers(state, endpoint_scaling_model(relperm), Val(ph))
    if ph == :wog
        for c in ix
            @inbounds update_three_phase_relperm!(kr, relperm, phases, s, c, state.ConnateWater[c], scalers)
        end
    else
        if ph == :wo
            kr1, kr2 = relperm.krw, relperm.krow
        elseif ph == :og
            kr1, kr2 = relperm.krog, relperm.krg
        else
            @assert ph == :wg
            kr1, kr2 = relperm.krw, relperm.krg
        end
        for c in ix
            @inbounds update_two_phase_relperm!(kr, relperm, kr1, kr2, phases, s, c, scalers)
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

    if !isnothing(scalers)
        swcon = scalers.w[1, c]
    end
    Krw, Krow, Krog, Krg = get_three_phase_relperms(relperm, c, krw, krow, krog, krg, swcon, scalers)

    sw = s[w, c]
    so = s[o, c]
    sg = s[g, c]


    kr[w, c] = Krw(sw)
    kr[o, c] = three_phase_oil_relperm(Krow(so), Krog(so), swcon, sg, sw)
    kr[g, c] = Krg(sg)
end

Base.@propagate_inbounds @inline function update_two_phase_relperm!(kr, relperm, krw, krn, phase_ind, s, c, scalers)
    w, n = phase_ind
    reg = region(relperm.regions, c)
    krw = table_by_region(krw, reg)
    krn = table_by_region(krn, reg)

    Krw, Krow = get_two_phase_relperms(relperm, c, krw, krn, scalers)

    sw = s[w, c]
    sn = s[n, c]

    kr[w, c] = Krw(sw)
    kr[n, c] = Krow(sn)
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