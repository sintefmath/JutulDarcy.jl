using MultiComponentFlash
export MultiPhaseCompositionalSystemLV
export StandardVolumeSource, VolumeSource, MassSource

const MINIMUM_COMPOSITIONAL_SATURATION = 1e-10

@inline function is_pure_single_phase(s_immiscible)
    return s_immiscible > 1.0 - MINIMUM_COMPOSITIONAL_SATURATION
end

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

function convergence_criterion(model::CompositionalModel, storage, eq::ConservationLaw, eq_s, r; dt = 1.0, update_report = missing)
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


function convergence_criterion(model::SimulationModel{<:Any, S}, storage, eq::ConservationLaw, eq_s, r; dt = 1.0, update_report = missing) where S<:MultiPhaseCompositionalSystemLV
    sys = model.system
    state = storage.state
    active = active_entities(model.domain, Cells())
    nc = number_of_components(sys)
    get_sat(ph) = as_value(view(state.Saturations, ph, :))
    get_density(ph) = as_value(view(state.PhaseMassDensities, ph, :))

    dz = compositional_increment(model, state, update_report)
    dp_abs, dp_rel = pressure_increments(model, state, update_report)

    dz0 = dz
    if has_other_phase(sys)
        a, l, v = phase_indices(sys)
        sw = get_sat(a)
        water_density = get_density(a)
        dsw_max = immiscible_increment(model, state, update_report)
        dz = max(dz, dsw_max)
    else
        l, v = phase_indices(sys)
        sw = nothing
        water_density = nothing
    end
    liquid_density = get_density(l)
    vapor_density = get_density(v)

    sl = get_sat(l)
    sv = get_sat(v)
    vol = as_value(state.FluidVolume)

    w = map(x -> x.mw, sys.equation_of_state.mixture.properties)
    e = compositional_criterion(state, dt, active, r, nc, w, sl, liquid_density, sv, vapor_density, sw, water_density, vol)
    names = model.system.components
    R = (
        CNV = (errors = e, names = names),
        increment_dp_abs = (errors = (dp_abs, ), names = (raw"Δp", ), ),
        increment_dp_rel = (errors = (dp_rel, ), names = (raw"Δp", ), ),
        increment_dz = (errors = (dz, ), names = (raw"Δz", ), )
        )
    return R
end

function compositional_increment(model, state, update_report::Missing)
    return 1.0
end

function pressure_increments(model, state, update_report::Missing)
    max_p = maximum(value, state.Pressure)
    return (max_p, 1.0)
end

function pressure_increments(model, state, update_report)
    max_p = maximum(value, state.Pressure)
    dp = update_report[:Pressure]
    dp_abs = dp.max
    dp_rel = dp_abs/max_p
    return (dp_abs, dp_rel)
end

function compositional_increment(model, state, update_report)
    return update_report[:OverallMoleFractions].max_scaled
end

function immiscible_increment(model, state, ::Missing)
    return 1.0
end

function immiscible_increment(model, state, update_report)
    return update_report[:ImmiscibleSaturation].max
end

function compositional_residual_scale(cell, dt, w, sl, liquid_density, sv, vapor_density, sw, water_density, vol)
    if isnothing(sw)
        sw_i = 0.0
    else
        sw_i = sw[cell]
    end

    if sw_i > 1.0 - 1e-4
        scale = 0.0
    else
        # The convergence criterion is taken to be dimensionless in a similar
        # manner to the standard CNV type criterion. We scale everything by the
        # immiscible phase so that the convergence criterion is more relaxed in
        # cells that are close to or completely filled with the immiscible
        # phase.
        scale_lv = 1.0 - sw_i
        sl_scaled = sl[cell]/scale_lv
        denl = liquid_density[cell]
        denv = vapor_density[cell]
        total_density = denl * sl_scaled + denv * (1.0 - sl_scaled)
        scale = dt * mean(w) * (scale_lv / (vol[cell] * max(total_density, 1e-3)))
    end
    return scale
end

function compositional_criterion(state, dt, active, r, nc, w, sl, liquid_density, sv, vapor_density, sw, water_density, vol)
    e = fill(-Inf, nc)
    for (ix, i) in enumerate(active)
        scaling = compositional_residual_scale(i, dt, w, sl, liquid_density, sv, vapor_density, sw, water_density, vol)
        for c in 1:(nc-1)
            val = scaling*abs(r[c, ix])/w[c]
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

function compositional_criterion(state, dt, active, r, nc, w, sl, liquid_density, sv, vapor_density, sw::Nothing, water_density::Nothing, vol)
    e = zeros(nc)
    for (ix, i) in enumerate(active)
        scaling = compositional_residual_scale(i, dt, w, sl, liquid_density, sv, vapor_density, sw, water_density, vol)
        for c in 1:nc
            e[c] = max(e[c], scaling*abs(r[c, ix])/w[c])
        end
    end
    return e
end
