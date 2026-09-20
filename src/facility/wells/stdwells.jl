const SimpleWellDomain = DiscretizedDomain{<:SimpleWell}
const SimpleWellFlowModel = SimulationModel{<:SimpleWellDomain, <:MultiPhaseSystem}

const SimpleWellModel = SimulationModel{<:SimpleWellDomain, <:Any}

struct WellMassFractions <: FractionVariables end
Jutul.need_default_primary(model, ::WellMassFractions) = false
Jutul.absolute_increment_limit(::WellMassFractions) = 0.15

function Jutul.default_values(model, mf::WellMassFractions)
    nc = values_per_entity(model, mf)
    return [1/nc for _ in 1:nc]
end

struct SimpleWellSystem{T, P, Ref} <: MultiPhaseSystem
    ncomp::Int
    phases::P
    c::Float64
    rho_ref::T
end

const StandardWellFlowModel = SimulationModel{<:SimpleWellDomain, <:SimpleWellSystem}

function SimpleWellSystem(ncomp, phases; c = 1e-8,
        reference_densities = ones(ncomp),
        reference_phase_index = get_reference_phase_index(phases))
    reference_densities = tuple(reference_densities...)
    return SimpleWellSystem{typeof(reference_densities), typeof(phases), reference_phase_index}(
        ncomp, phases, c, reference_densities)
end

function SimpleWellSystem(system; reference_phase_index = get_reference_phase_index(system), kwarg...)
    rho = reference_densities(system)
    ncomp = number_of_components(system)
    phases = get_phases(system)
    return SimpleWellSystem(ncomp, phases;
        reference_densities = rho,
        reference_phase_index = reference_phase_index, kwarg...)
end

number_of_components(s::SimpleWellSystem) = s.ncomp
# number_of_phases(s::SimpleWellSystem) = s.ncomp
reference_densities(s::SimpleWellSystem) = s.rho_ref
reference_densities(s::Symbol) = reference_densities(Val(s))

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
    if model_is_thermal(model)
        prm[:WellIndicesThermal] = WellIndicesThermal()
    end
    prm[:PerforationGravityDifference] = PerforationGravityDifference()
end

function Jutul.initialize_extra_state_fields!(state, d::DiscretizedDomain, m::SimpleWellFlowModel; T = Float64)
    if well_has_explicit_pressure_drop(m)
        nc = count_entities(d, Perforations())
        state[:ConnectionPressureDrop] = zeros(T, nc)
    end
end

function Jutul.select_minimum_output_variables!(vars, domain::DiscretizedDomain, model::SimpleWellFlowModel)
    push!(vars, :PhaseMassDensities)
    push!(vars, :Saturations)
    return vars
end

@inline get_reference_phase_index(::SimpleWellSystem{T, P, Ref}) where {T, P, Ref} = Ref

well_has_explicit_pressure_drop(m::SimpleWellFlowModel) = well_has_explicit_pressure_drop(physical_representation(m.domain))
well_has_explicit_pressure_drop(w::SimpleWell) = w.explicit_dp

function connection_pressure_drop_mask_values(mask)
    return mask.values
end

function connection_pressure_drop_mask_values(::Nothing)
    return nothing
end

function connection_pressure_drop_mask_values(mask, context)
    return Adapt.adapt(context, mask.values)
end

function connection_pressure_drop_mask_values(::Nothing, context)
    return nothing
end

function update_before_step_well!(well_state,
        well_model::SimpleWellFlowModel,
        res_state, res_model, ctrl, mask;
        update_explicit = true,
        backend_well_state = well_state,
        backend_well_model = well_model,
        backend_reservoir_state = res_state,
        backend_reservoir_model = res_model)
    if well_has_explicit_pressure_drop(well_model) && update_explicit
        if res_model.context isa Jutul.KernelAbstractionsContext
            update_connection_pressure_drop_backend!(
                well_state.ConnectionPressureDrop,
                backend_well_state, backend_well_model,
                backend_reservoir_state, backend_reservoir_model,
                ctrl, mask)
        else
            update_connection_pressure_drop!(
                well_state.ConnectionPressureDrop,
                well_state, well_model, res_state, res_model, ctrl, mask)
        end
    end
    return nothing
end

function update_connection_pressure_drop!(
        dp, well_state, well_model, res_state, res_model,
        ctrl::InjectorControl, mask)
    phases = ctrl.phases
    perf = physical_representation(well_model.domain).perforations
    res_cells = perf.reservoir
    gdz = well_state.PerforationGravityDifference
    ρ = res_state.PhaseMassDensities
    return update_injector_connection_pressure_drop!(
        dp, phases, res_cells, gdz, ρ)
end

function update_injector_connection_pressure_drop!(
        dp, phases, res_cells, gdz, ρ)
    # Traverse down the well, using the phase notion encoded in ctrl and then
    # just accumulate pressure drop as we go assuming no cross flow.
    T = eltype(dp)
    dp_current = zero(T)
    gdz_current = zero(T)
    for i in eachindex(dp)
        rc = res_cells[i]
        gdz_next = value(gdz[i])

        Δgdz = gdz_next - gdz_current
        # Mixture density along well bore
        local_density = zero(T)
        for (ph, mix) in phases
            local_density += convert(T, mix)*value(ρ[ph, rc])
        end
        dp_current += local_density*Δgdz
        dp[i] = value(dp_current)

        # Onto next one
        gdz_current = gdz_next
    end
    return dp
end

function update_connection_pressure_drop!(
        dp, well_state, well_model, res_state, res_model, ctrl, mask)
    # Well is either disabled or producing. Loop over well from the bottom,
    # aggregating mixture density as we go. Then traverse down from the top and
    # accumulate the actual pressure drop due to hydrostatic assumptions.
    well = physical_representation(well_model)
    perf = well.perforations
    res_cells = perf.reservoir
    gdz = well_state.PerforationGravityDifference
    WI = well_state.WellIndices
    ρ = res_state.PhaseMassDensities
    mob = res_state.PhaseMobilities
    p = res_state.Pressure
    well_pressure = well_state.Pressure
    mask_values = connection_pressure_drop_mask_values(mask)
    return update_producer_connection_pressure_drop!(
        dp, res_cells, gdz, WI, ρ, mob, p, well_pressure, mask_values)
end

function update_producer_connection_pressure_drop!(
        dp, res_cells, gdz, WI, ρ, mob, p, well_pressure, mask_values)
    T = eltype(dp)
    bhp = value(well_pressure[1])
    # Integrate up, adding weighted density into well bore and keeping track of
    # current weight. Initialize by a simple average first.
    mobility_density_sum = zero(T)
    mobility_sum = zero(T)
    for i in eachindex(dp)
        rc = res_cells[i]
        mobility_density_sum = zero(T)
        WI_i = value(WI[i])
        if !isnothing(mask_values)
            WI_i *= value(mask_values[i])
        end
        for ph in axes(ρ, 1)
            mob_ph = WI_i*value(mob[ph, rc])
            mobility_sum += mob_ph
            mobility_density_sum += value(ρ[ph, rc])*mob_ph
        end
        if i == 1
            dz = value(gdz[i])
            dp_prev = zero(T)
        else
            dz = value(gdz[i]) - value(gdz[i-1])
            dp_prev = dp[i-1]
        end
        est_density = mobility_density_sum/max(mobility_sum, convert(T, 1e-3))
        dp[i] = dp_prev + dz*est_density
    end
    current_density = current_weight = zero(T)
    for i in reverse(eachindex(dp))
        rc = res_cells[i]
        wi = value(WI[i])
        if !isnothing(mask_values)
            wi *= value(mask_values[i])
        end
        pot = abs(bhp + dp[i] - value(p[rc]))
        q_perf = wi*pot
        # Mixture density along well bore
        local_density = zero(T)
        local_weight = zero(T)
        for ph in axes(ρ, 1)
            λ = value(mob[ph, rc])
            weight_ph = q_perf*λ
            local_weight += weight_ph
            local_density += weight_ph*value(ρ[ph, rc])
        end
        current_weight += local_weight
        current_density += local_density
        if abs(current_weight) > zero(T)
            next_dp = current_density/current_weight
        else
            next_dp = zero(T)
        end
        # NOTE: This is a real hack - should be fixed higher up in the
        # initialization of the parameter.
        dp[i] = value(next_dp)
    end
    # Integrate down, using the mixture densities (temporarily stored in dp) to
    # calculate the pressure drop from the top.
    dp_current = zero(T)
    gdz_current = zero(T)
    for i in eachindex(dp)
        local_density = dp[i]
        gdz_next = value(gdz[i])
        Δgdz = gdz_next - gdz_current

        dp_current += local_density*Δgdz
        gdz_current = gdz_next

        dp[i] = value(dp_current)
    end
    return dp
end

function update_connection_pressure_drop_backend!(
        host_dp, well_state, well_model, res_state, res_model,
        ctrl::InjectorControl, mask)
    context = res_model.context
    dp = well_state.ConnectionPressureDrop
    perf = physical_representation(well_model.domain).perforations
    res_cells = perf.reservoir
    gdz = well_state.PerforationGravityDifference
    ρ = res_state.PhaseMassDensities
    phases = Tuple(ctrl.phases)
    function update_pressure_drop(_)
        update_injector_connection_pressure_drop!(
            dp, phases, res_cells, gdz, ρ)
        return nothing
    end
    Jutul.threaded_loop(update_pressure_drop, 1, context)
    Jutul.backend_copyto!(host_dp, dp)
    Jutul.synchronize(context)
    return nothing
end

function update_connection_pressure_drop_backend!(
        host_dp, well_state, well_model, res_state, res_model,
        ctrl, mask)
    context = res_model.context
    dp = well_state.ConnectionPressureDrop
    perf = physical_representation(well_model).perforations
    res_cells = perf.reservoir
    gdz = well_state.PerforationGravityDifference
    WI = well_state.WellIndices
    ρ = res_state.PhaseMassDensities
    mob = res_state.PhaseMobilities
    p = res_state.Pressure
    well_pressure = well_state.Pressure
    mask_values = connection_pressure_drop_mask_values(mask, context)
    function update_pressure_drop(_)
        update_producer_connection_pressure_drop!(
            dp, res_cells, gdz, WI, ρ, mob, p,
            well_pressure, mask_values)
        return nothing
    end
    Jutul.threaded_loop(update_pressure_drop, 1, context)
    Jutul.backend_copyto!(host_dp, dp)
    Jutul.synchronize(context)
    return nothing
end
