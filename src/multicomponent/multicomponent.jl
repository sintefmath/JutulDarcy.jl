using MultiComponentFlash
export MultiPhaseCompositionalSystemLV
export StandardVolumeSource, VolumeSource, MassSource

const MINIMUM_COMPOSITIONAL_SATURATION = 1e-10

include("variables/variables.jl")
include("utils.jl")
include("flux.jl")
include("sources.jl")
include("wells.jl")

function select_primary_variables!(S, system::CompositionalSystem, model)
    S[:Pressure] = Pressure(max_rel = 0.2, minimum = 1e5)
    S[:OverallMoleFractions] = OverallMoleFractions(dz_max = 0.1)
    if has_other_phase(system)
        S[:ImmiscibleSaturation] = ImmiscibleSaturation(ds_max = 0.2)
    end
end

function select_secondary_variables!(S, system::CompositionalSystem, model)
    select_default_darcy_secondary_variables!(S, model.domain, system, model.formulation)
    if has_other_phase(system)
        water_pvt = ConstMuBTable(101325.0, 1.0, 1e-18, 1e-3, 1e-20)
        set_secondary_variables!(model, PhaseViscosities = ThreePhaseLBCViscositiesLV(water_pvt),
                                        PhaseMassDensities = ThreePhaseCompositionalDensitiesLV(water_pvt))
    else
        set_secondary_variables!(model, PhaseViscosities = LBCViscosities(),
                                        PhaseMassDensities = TwoPhaseCompositionalDensities())
    end
    S[:LiquidMassFractions] = PhaseMassFractions(:liquid)
    S[:VaporMassFractions] = PhaseMassFractions(:vapor)
    S[:FlashResults] = FlashResults(system)
    S[:Saturations] = Saturations()
end

function select_parameters!(prm, system::CompositionalSystem, model)
    select_default_darcy_parameters!(prm, model.domain, system, model.formulation)
    prm[:Temperature] = Temperature()
end

function convergence_criterion(model::CompositionalModel, storage, eq::ConservationLaw, eq_s, r; dt = 1)
    tm = storage.state0.TotalMasses
    a = active_entities(model.domain, Cells())
    function scale(i)
        @inbounds c = a[i]
        t = 0.0
        @inbounds for i = 1:size(tm, 1)
            t += tm[i, c]
        end
        return t
    end
    @tullio max e[j] := abs(r[j, i]) * dt / scale(i)
    names = model.system.components
    R = Dict("CNV" => (errors = e, names = names))
    return R
end


function convergence_criterion(model::SimulationModel{<:Any, S}, storage, eq::ConservationLaw, eq_s, r; dt = 1) where S<:MultiPhaseCompositionalSystemLV
    tm = storage.state0.TotalMasses

    sys = model.system
    a = active_entities(model.domain, Cells())
    nc = number_of_components(sys)
    has_water = has_other_phase(sys)
    if has_water
        water = storage.state.ImmiscibleSaturation
    else
        water = nothing
    end

    e = compositional_criterion(dt, tm, a, r, nc, water)
    names = model.system.components
    R = Dict("CNV" => (errors = e, names = names))
    return R
end


function compositional_residual_scale(cell, dt, tm)
    t = 0.0
    @inbounds for i = 1:size(tm, 1)
        t += tm[i, cell]
    end
    return dt/t
end

function compositional_criterion(dt, total_mass0, active, r, nc, water)
    e = fill(-Inf, nc)
    worst = Vector{Int64}(undef, nc)
    for (ix, i) in enumerate(active)
        s = compositional_residual_scale(i, dt, total_mass0)
        sw = value(water[i])
        sc = (1 - sw)*s
        for c in 1:(nc-1)
            val =  sc*abs(r[c, ix])
            if val > e[c]
                e[c] = val
                worst[c] = i
            end
        end
        valw = s*abs(r[end, ix])
        if valw > e[end]
            e[end] = valw
            worst[end] = i
        end
    end
    # @debug "Worst cells" e worst
    return e
end

function compositional_criterion(dt, total_mass0, active, r, nc, water::Nothing)
    e = zeros(nc)
    for (ix, i) in enumerate(active)
        s = compositional_residual_scale(i, dt, total_mass0)
        for c in 1:nc
            e[c] = max(e[c], s*abs(r[c, ix]))
        end
    end
    return e
end
