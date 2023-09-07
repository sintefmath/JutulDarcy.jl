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
    S[:Pressure] = Pressure()
    S[:OverallMoleFractions] = OverallMoleFractions(dz_max = 0.1)
    if has_other_phase(system)
        S[:ImmiscibleSaturation] = ImmiscibleSaturation(ds_max = 0.2)
    end
end

function select_secondary_variables!(S, system::CompositionalSystem, model)
    select_default_darcy_secondary_variables!(S, model.domain, system, model.formulation)
    if has_other_phase(system)
        water_pvt = ConstMuBTable(DEFAULT_MINIMUM_PRESSURE, 1.0, 1e-18, 1e-3, 1e-20)
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
        @inbounds for i in axes(tm, 1)
            t += tm[i, c]
        end
        return t
    end
    @tullio max e[j] := abs(r[j, i]) * dt / scale(i)
    names = model.system.components
    R = (CNV = (errors = e, names = names), )
    return R
end


function convergence_criterion(model::SimulationModel{<:Any, S}, storage, eq::ConservationLaw, eq_s, r; dt = 1) where S<:MultiPhaseCompositionalSystemLV
    tm = storage.state0.TotalMasses

    sys = model.system
    a = active_entities(model.domain, Cells())
    nc = number_of_components(sys)
    has_water = has_other_phase(sys)
    if has_water
        a, l, v = phase_indices(sys)
        water = as_value(storage.state.ImmiscibleSaturation)
        water_density = as_value(view(storage.state.PhaseMassDensities, a, :))
    else
        water = nothing
        water_density = nothing
    end
    vol = as_value(storage.state.FluidVolume)

    w = map(x -> x.mw, sys.equation_of_state.mixture.properties)
    e = compositional_criterion(dt, tm, a, r, nc, w, water, water_density, vol)
    names = model.system.components
    R = (CNV = (errors = e, names = names), )
    return R
end


function compositional_residual_scale(cell, dt, w, tm)
    t = 0.0
    @inbounds for i = axes(tm, 1)
        t += tm[i, cell]/w[i]
    end
    return dt/t
end


function compositional_criterion(dt, total_mass0, active, r, nc, w, water, water_density, vol)
    e = fill(-Inf, nc)
    for (ix, i) in enumerate(active)
        s = compositional_residual_scale(i, dt, w, total_mass0)
        sw = value(water[i])
        sc = (1 - sw)*s
        for c in 1:(nc-1)
            val = sc*abs(r[c, ix])/w[c]
            if val > e[c]
                e[c] = val
            end
        end
        valw = dt*abs(r[end, ix])/(water_density[i]*vol[i])
        if valw > e[end]
            e[end] = valw
        end
    end
    return e
end

function compositional_criterion(dt, total_mass0, active, r, nc, w, water::Nothing, ::Nothing, vol)
    e = zeros(nc)
    for (ix, i) in enumerate(active)
        s = compositional_residual_scale(i, dt, w, total_mass0)
        for c in 1:nc
            e[c] = max(e[c], s*abs(r[c, ix])/w[c])
        end
    end
    return e
end
