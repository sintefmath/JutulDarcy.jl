const SimpleWellDomain = DiscretizedDomain{<:SimpleWell}
const SimpleWellFlowModel = SimulationModel{<:SimpleWellDomain, <:MultiPhaseSystem}

struct WellMassFractions <: FractionVariables end
Jutul.need_default_primary(model, ::WellMassFractions) = false
function Jutul.default_values(model, mf::WellMassFractions)
    nc = values_per_entity(model, mf)
    return [1/nc for _ in 1:nc]
end

struct SimpleWellEquation <: JutulEquation end

Jutul.local_discretization(::SimpleWellEquation, i) = nothing

function Jutul.number_of_equations_per_entity(model::SimpleWellFlowModel, ::SimpleWellEquation)
    sys = model.system
    return number_of_components(sys)
end

function values_per_entity(model, v::WellMassFractions)
    sys = model.system
    return number_of_components(sys)
end

function select_primary_variables!(model::SimpleWellFlowModel)
    pvars = model.primary_variables
    pvars[:Pressure] = Pressure()
    pvars[:MassFractions] = WellMassFractions()
end

function select_secondary_variables!(model::SimpleWellFlowModel)

end

function select_equations!(model::SimpleWellFlowModel)
    model.equations[:mass_conservation] = SimpleWellEquation()
end

function select_parameters!(model::SimpleWellFlowModel)
    prm = model.parameters
    prm[:WellIndices] = WellIndices()
    prm[:PerforationGravityDifference] = PerforationGravityDifference()
end

function select_minimum_output_variables!(model::SimpleWellFlowModel)
    outputs = model.output_variables
    for k in keys(model.primary_variables)
        push!(outputs, k)
    end
    return model
end

function Jutul.initialize_extra_state_fields!(state, d::DiscretizedDomain, m::SimpleWellFlowModel)
    nc = count_entities(d, Perforations())
    state[:ConnectionPressureDrop] = zeros(nc)
end

well_has_explicit_pressure_drop(m::SimpleWellFlowModel) = well_has_explicit_pressure_drop(m.domain.grid)
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
    perf = well_model.domain.grid.perforations
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
    perf = well_model.domain.grid.perforations
    res_cells = perf.reservoir
    WI = perf.WI

    gdz = perf.gdz

    ρ = as_value(res_state.PhaseMassDensities)
    kr = as_value(res_state.RelativePermeabilities)
    mu = as_value(res_state.PhaseViscosities)

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
            λ = kr[ph, rc]/mu[ph, rc]
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
