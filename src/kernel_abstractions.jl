# Structural adaptation for the compact topology used by the single-model
# reservoir discretization. More complex property types can add equivalent
# Adapt rules without changing the execution context.
function Adapt.adapt_structure(to, g::MinimalTPFATopology)
    return MinimalTPFATopology(g.nc, Adapt.adapt(to, g.neighborship))
end

function _ka_cnv_mb_errors(r, pore_volume, density, dt, ::Val{N}, ctx) where N
    T = eltype(r)
    cnv = Jutul.backend_allocate(ctx, T, N)
    mb = Jutul.backend_allocate(ctx, T, N)
    nc = length(pore_volume)
    function reduce_phase(phase)
        local_cnv = zero(T)
        local_mb = zero(T)
        density_sum = zero(T)
        total_pore_volume = zero(T)
        for cell in 1:nc
            @inbounds begin
                pv = pore_volume[cell]
                rho = value(density[phase, cell])
                residual = r[phase, cell]
                total_pore_volume += pv
                local_mb += residual
                density_sum += abs(rho)
                local_cnv = max(local_cnv, dt*abs(residual)/(rho*pv))
            end
        end
        average_density = density_sum/nc
        @inbounds begin
            cnv[phase] = local_cnv
            mb[phase] = (dt/total_pore_volume)*abs(local_mb)/average_density
        end
    end
    Jutul.threaded_loop(reduce_phase, N, ctx)
    return (Tuple(Jutul.backend_to_host(cnv)), Tuple(Jutul.backend_to_host(mb)))
end

function convergence_criterion(
        model::SimulationModel{D, S, F, C}, storage,
        eq::ConservationLaw{:TotalMasses}, eq_s, r;
        dt = 1, update_report = missing
    ) where {D, S<:MultiPhaseSystem, F<:JutulFormulation,
             C<:Jutul.KernelAbstractionsContext}
    map = Jutul.global_map(model.domain)
    active(x) = Jutul.active_view(x, map, for_variables = false)
    pore_volume = active(storage.state.FluidVolume)
    density = active(storage.state.PhaseMassDensities)
    nph = number_of_phases(model.system)
    cnv, mb = _ka_cnv_mb_errors(r, pore_volume, density, dt, Val(nph), model.context)
    dp_abs, dp_rel = pressure_increments(model, storage.state, update_report)
    if ismissing(update_report)
        ds_max = 1.0
    elseif haskey(update_report, :Saturations)
        ds_max = update_report[:Saturations].max
    else
        ds_max = 0.0
    end
    names = phase_names(model.system)
    return (
        CNV = (errors = cnv, names = names),
        MB = (errors = mb, names = names),
        increment_dp_abs = (errors = (dp_abs/1e6, ), names = (raw"Δp (abs, MPa)", ), ),
        increment_dp_rel = (errors = (dp_rel, ), names = (raw"Δp (rel)", ), ),
        increment_saturation = (errors = (ds_max, ), names = (raw"Δs", ), )
    )
end

function pressure_increments(model::SimulationModel{D, S, F, C}, state,
        update_report::Missing) where {D, S, F<:JutulFormulation,
                                       C<:Jutul.KernelAbstractionsContext}
    max_pressure = Jutul.backend_maximum_value(model.context, state.Pressure)
    return (max_pressure, 1.0)
end

function pressure_increments(model::SimulationModel{D, S, F, C}, state,
        update_report) where {D, S, F<:JutulFormulation,
                              C<:Jutul.KernelAbstractionsContext}
    max_pressure = Jutul.backend_maximum_value(model.context, state.Pressure)
    dp_abs = update_report[:Pressure].max
    return (dp_abs, dp_abs/max_pressure)
end

# Forces stay as small CPU tuples/vectors. Launching one kernel per force avoids
# races when two entries address the same cell and avoids adapting a dynamic
# Vector container into a model kernel.
function Jutul.apply_forces_to_equation!(acc, storage,
        model::SimulationModel{D, S, F, C}, eq::ConservationLaw,
        eq_s, sources::V, time
    ) where {D, S<:MultiPhaseSystem, F<:JutulFormulation,
             C<:Jutul.KernelAbstractionsContext,
             I, Q, FF, V<:AbstractVector{SourceTerm{I, Q, FF}}}
    state = storage.state
    kr = haskey(state, :RelativePermeabilities) ? state.RelativePermeabilities : 1.0
    viscosity = state.PhaseViscosities
    rho_surface = reference_densities(model.system)
    nph = size(acc, 1)
    map = Jutul.global_map(model.domain)
    for source in sources
        function apply_source(_)
            cell = Jutul.full_cell(source.cell, map)
            for phase in 1:nph
                q = phase_source(cell, source, rho_surface[phase], kr, viscosity, phase)
                @inbounds acc[phase, source.cell] -= q
            end
        end
        Jutul.threaded_loop(apply_source, 1, model.context)
    end
    return acc
end

function Jutul.apply_forces_to_equation!(acc, storage,
        model::SimulationModel{D, S, F, C}, eq::ConservationLaw{:TotalMasses},
        eq_s, conditions::V, time
    ) where {D, S<:MultiPhaseSystem, F<:JutulFormulation,
             C<:Jutul.KernelAbstractionsContext,
             V<:AbstractVector{<:FlowBoundaryCondition}}
    state = storage.state
    system = model.system
    map = Jutul.global_map(model.domain)
    for condition in conditions
        function apply_condition(_)
            flux = compute_bc_mass_fluxes(system, condition, map, state)
            cell = condition.cell
            for phase in axes(acc, 1)
                @inbounds acc[phase, cell] += flux[phase]
            end
        end
        Jutul.threaded_loop(apply_condition, 1, model.context)
    end
    return acc
end
