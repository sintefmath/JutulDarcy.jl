const SimpleWellDomain = DiscretizedDomain{<:SimpleWell}
const SimpleWellFlowModel = SimulationModel{<:SimpleWellDomain, <:MultiPhaseSystem}

struct WellMassFractions <: FractionVariables end
Jutul.need_default_primary(model, ::WellMassFractions) = false
Jutul.absolute_increment_limit(::WellMassFractions) = 0.15

function Jutul.default_values(model, mf::WellMassFractions)
    nc = values_per_entity(model, mf)
    return [1/nc for _ in 1:nc]
end

struct SimpleWellSystem{T, P} <: MultiPhaseSystem
    ncomp::Int
    phases::P
    c::Float64
    rho_ref::T
end

const StandardWellFlowModel = SimulationModel{<:SimpleWellDomain, <:SimpleWellSystem}

function SimpleWellSystem(ncomp, phases; c = 1e-8, reference_densities = ones(ncomp))
    reference_densities = tuple(reference_densities...)
    return SimpleWellSystem(ncomp, phases, c, reference_densities)
end

function SimpleWellSystem(system; kwarg...)
    rho = reference_densities(system)
    ncomp = number_of_components(system)
    phases = get_phases(system)
    return SimpleWellSystem(ncomp, phases; reference_densities = rho, kwarg...)
end

number_of_components(s::SimpleWellSystem) = s.ncomp
# number_of_phases(s::SimpleWellSystem) = s.ncomp
reference_densities(s::SimpleWellSystem) = s.rho_ref

function flash_wellstream_at_surface(var, well_model, system::SimpleWellSystem, well_state, rhoS, cond = default_surface_cond())
    X = well_state.MassFractions
    vol = X./rhoS
    volfrac = vol./sum(vol)
    return (rhoS, volfrac)
end

function Jutul.values_per_entity(model, v::WellMassFractions)
    sys = model.system
    return number_of_components(sys)
end

function Jutul.select_primary_variables!(pvars, s::SimpleWellSystem, model::SimpleWellFlowModel)
    pvars[:Pressure] = Pressure(max_rel = Inf)
    pvars[:MassFractions] = WellMassFractions()
end

function Jutul.select_secondary_variables!(S, system::SimpleWellSystem, model::SimpleWellFlowModel)
    S[:TotalMasses] = TotalMasses()
end

function select_parameters!(prm, s::SimpleWellDomain, model::SimpleWellFlowModel)
    prm[:FluidVolume] = FluidVolume()
    prm[:WellIndices] = WellIndices()
    prm[:PerforationGravityDifference] = PerforationGravityDifference()
end

function Jutul.initialize_extra_state_fields!(state, d::DiscretizedDomain, m::SimpleWellFlowModel)
    if well_has_explicit_pressure_drop(m)
        nc = count_entities(d, Perforations())
        state[:ConnectionPressureDrop] = zeros(nc)
    end
end

well_has_explicit_pressure_drop(m::SimpleWellFlowModel) = well_has_explicit_pressure_drop(physical_representation(m.domain))
well_has_explicit_pressure_drop(w::SimpleWell) = w.explicit_dp

function update_before_step_well!(well_state, well_model::SimpleWellFlowModel, res_state, res_model, ctrl)
    if well_has_explicit_pressure_drop(well_model)
        dp = well_state.ConnectionPressureDrop
        update_connection_pressure_drop!(dp, well_state, well_model, res_state, res_model, ctrl)
    end
end

function update_connection_pressure_drop!(dp, well_state, well_model, res_state, res_model, ctrl::InjectorControl)
    # Traverse down the well, using the phase notion encoded in ctrl and then
    # just accumulate pressure drop as we go assuming no cross flow
    phases = ctrl.phases
    perf = physical_representation(well_model.domain).perforations
    res_cells = perf.reservoir
    gdz = perf.gdz

    ρ = as_value(res_state.PhaseMassDensities)

    dp_current = 0.0
    gdz_current = 0.0
    for i in eachindex(dp)
        rc = res_cells[i]
        gdz_next = gdz[i]

        Δgdz = gdz_next - gdz_current
        # Mixture density along well bore
        local_density = 0
        for (ph, mix) in phases
            local_density += mix*ρ[ph, rc]
        end
        dp_current += local_density*Δgdz
        dp[i] = dp_current

        # Onto next one
        gdz_current = gdz_next
    end
end

function update_connection_pressure_drop!(dp, well_state, well_model, res_state, res_model, ctrl)
    # Well is either disabled or producing. Loop over well from the bottom,
    # aggregating mixture density as we go. Then traverse down from the top and
    # accumulate the actual pressure drop due to hydrostatic assumptions.
    perf = physical_representation(well_model).perforations
    res_cells = perf.reservoir
    WI = perf.WI

    gdz = perf.gdz

    ρ = as_value(res_state.PhaseMassDensities)
    mob = as_value(res_state.PhaseMobilities)
    # kr = as_value(res_state.RelativePermeabilities)
    # mu = as_value(res_state.PhaseViscosities)

    # Integrate up, adding weighted density into well bore and keeping track of
    # current weight
    current_weight = 0.0
    current_density = 0.0
    for i in reverse(eachindex(dp))
        rc = res_cells[i]
        wi = WI[i]

        # Mixture density along well bore
        local_density = 0
        local_weight = 0
        for ph in axes(ρ, 1)
            λ = mob[ph, rc]
            weight_ph = wi*λ
            local_weight += weight_ph
            local_density += weight_ph*ρ[ph, rc]
        end
        current_weight += local_weight
        current_density += local_density
        dp[i] = current_density/current_weight
    end
    # Integrate down, using the mixture densities (temporarily stored in dp) to
    # calculate the pressure drop from the top.
    dp_current = 0.0
    gdz_current = 0.0
    for i in eachindex(dp)
        local_density = dp[i]
        gdz_next = gdz[i]
        Δgdz = gdz_next - gdz_current

        dp_current += local_density*Δgdz
        gdz_current += gdz_next

        dp[i] = dp_current
    end
end
