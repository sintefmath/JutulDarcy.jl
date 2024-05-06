function state_pair(storage_g, storage_l, current = true)
    # storage_g = simulator.simulator.storage
    # storage_l = sim.storage
    if current
        return (storage_g.state, storage_l.state)
    else
        return (storage_g.state0, storage_l.state0)
    end
end

function update_subdomain_from_global(simulator, sim, i; kwarg...)
    # Add a nesting so that we can dispatch on the model type
    update_subdomain_from_global(sim.model, simulator, sim, i; kwarg...)
end

function update_subdomain_from_global(inner_model, simulator, sim, i; current = true, kwarg...)
    s_l = sim.storage
    s_g = simulator.storage
    m = sim.model
    Ω = m.domain
    M = global_map(Ω)

    (state_g, state_l) = state_pair(s_g, s_l, current)
    var_def_g = s_g.variable_definitions
    var_def_l = s_l.variable_definitions

    @tic "primary to local" primary_to_local!(state_l, state_g, var_def_g, m, M)
    handle_secondary(state_l, state_g, var_def_l, var_def_g, m, M; kwarg...)
end

function update_subdomain_from_global(inner_model::MultiModel, simulator, sim, i; current = true, kwarg...)
    storage_g = simulator.storage
    storage_l = sim.storage
    models = inner_model.models

    for k in keys(models)
        m = models[k]
        s_g, s_l = storage_g[k], storage_l[k]
        var_def_g = s_g.variable_definitions
        var_def_l = s_l.variable_definitions
        (state_g, state_l) = state_pair(s_g, s_l, current)
        Ω = m.domain
        M = global_map(Ω)
        @tic "primary to local" primary_to_local!(state_l, state_g, var_def_g, m, M)
        handle_secondary(state_l, state_g, var_def_l, var_def_g, m, M; kwarg...)
    end
end

function handle_secondary(state_l, state_g, defs_l, defs_g, m, M; transfer_secondary = false)
    if transfer_secondary
        @tic "secondary to local" secondary_to_local!(state_l, state_g, defs_g, m, M)
    else
        @tic "secondary variables" Jutul.update_secondary_variables_state!(state_l, m, defs_l.secondary_variables)
    end
end

function update_global_from_subdomain(simulator, sim, i; kwarg...)
    # Add a nesting so that we can dispatch on the model type
    update_global_from_subdomain(sim.model, simulator, sim, i; kwarg...)
end

function update_global_from_subdomain(inner_model, simulator, sim, i; secondary = true)
    model = sim.model
    Ω = model.domain
    M = Ω.global_map
    (state_g, state_l) = state_pair(simulator.storage, sim.storage, true)
    defs = sim.storage.variable_definitions
    @tic "primary to global" primary_to_global!(state_g, state_l, defs, model, M)
    if secondary
        @tic "secondary to global" secondary_to_global!(state_g, state_l, defs, model, M)
    end
end

function update_global_from_subdomain(inner_model::MultiModel, simulator, sim, i; secondary = true)
    storage_g = simulator.storage
    storage_l = sim.storage
    models = inner_model.models
    for k in keys(models)
        m = models[k]
        M = global_map(m.domain)
        s_g, s_l = storage_g[k], storage_l[k]
        defs = s_g.variable_definitions

        (state_g, state_l) = state_pair(s_g, s_l, true)
        @tic "primary to global" primary_to_global!(state_g, state_l, defs, m, M)
        if secondary
            @tic "secondary to global" secondary_to_global!(state_g, state_l, defs, m, M)
        end
    end
end

function variables_to_global!(state_global, state, M, variables)
    for k in variables
        v_global = state_global[k]
        v_local = state[k]
        transfer_to_global!(v_global, v_local, M)
    end
end

function variables_to_local!(state, state_global, M, variables)
    for k in variables
        v_global = state_global[k]
        v_local = state[k]
        transfer_to_local!(v_local, v_global, M)
    end
end

function primary_to_global!(state_global, state, defs, model, M)
    variables_to_global!(state_global, state, M, keys(defs.primary_variables))
    cfg = :WellGroupConfiguration
    if haskey(state_global, cfg)
        transfer_to_global!(state_global[cfg], state[cfg], M)
    end
    dp = :ConnectionPressureDrop
    if haskey(state_global, dp)
        # TODO: This part should not be needed.
        # transfer_to_local!(state_global[dp], state[dp], M)
    end
end

function secondary_to_global!(state_global, state, defs, model, M)
    variables_to_global!(state_global, state, M, keys(defs.secondary_variables))
end

function secondary_to_local!(state_global, state, defs, model, M)
    variables_to_local!(state_global, state, M, keys(defs.secondary_variables))
end

function primary_to_local!(state_global, state, defs, model, M)
    variables_to_local!(state_global, state, M, keys(defs.primary_variables))
    cfg = :WellGroupConfiguration
    if haskey(state_global, cfg)
        transfer_to_local!(state_global[cfg], state[cfg], M)
    end
    dp = :ConnectionPressureDrop
    if haskey(state_global, dp)
        transfer_to_local!(state_global[dp], state[dp], M)
    end
end

function previous_to_local!(state_global, state, defs, model, M)
    # variables_to_local!(state_global, state, M, model.output_variables)
    variables_to_local!(state_global, state, M, keys(defs.primary_variables))
    variables_to_local!(state_global, state, M, keys(defs.secondary_variables))
end

"Transfer to global - trivial map"
function transfer_to_global!(glob, loc, M)
    for i in eachindex(loc)
        @inbounds glob[i] = loc[i]
    end
end

function transfer_to_global!(glob::W, loc::W, M) where W<:JutulDarcy.WellGroupConfiguration
    transfer_well_config!(glob, loc, from_local = true)
end

function transfer_to_local!(loc::W, glob::W, M) where W<:JutulDarcy.WellGroupConfiguration
    # @info "Transferring wells to local"
    transfer_well_config!(loc, glob, from_local = false)
end

function transfer_well_config!(loc::W, glob::W; from_local::Bool) where W<:JutulDarcy.WellGroupConfiguration
    for k in keys(loc.operating_controls)
        if haskey(glob.operating_controls, k)
            loc_ctrl = loc.operating_controls[k]
            # skip_transfer = !(from_local && loc_ctrl isa DisabledControl)
            skip_transfer = false
            if skip_transfer
                # Don't propagate disabled controls over. Let the outer solver decide.
            else
                loc.operating_controls[k] = glob.operating_controls[k]
            end
        end
    end
    for k in keys(loc.requested_controls)
        if haskey(glob.requested_controls, k)
            loc.requested_controls[k] = glob.requested_controls[k]
        end
    end
    for k in keys(loc.limits)
        if haskey(glob.limits, k)
            loc.limits[k] = glob.limits[k]
        end
    end
end

"Transfer to global - vector, non-trivial map"
function transfer_to_global!(glob::AbstractVector, loc::AbstractVector, M::Jutul.FiniteVolumeGlobalMap{R}) where R
    @inbounds for c in eachindex(loc)
        if !Jutul.cell_is_boundary(c, M)
            @inbounds gc = Jutul.global_cell(c, M)::R
            @inbounds glob[gc] = loc[c]
        end
    end
end

"Transfer to global - matrix, non-trivial map"
function transfer_to_global!(glob::AbstractMatrix, loc::AbstractMatrix, M::Jutul.FiniteVolumeGlobalMap{R}) where R
    @inbounds for c in axes(loc, 2)
        if !Jutul.cell_is_boundary(c, M)
            @inbounds gc = Jutul.global_cell(c, M)::R
            @inbounds for d in axes(loc, 1)
                val = loc[d, c]
                glob[d, gc] = val
            end
        end
    end
end

"Transfer - trivial mapping"
function transfer_to_local!(loc, glob, M)
    for i in eachindex(loc)
        @inbounds loc[i] = glob[i]
    end
end

function transfer_to_local!(loc::Dict, glob::Dict, M::Jutul.TrivialGlobalMap)
    for k in keys(loc)
        loc[k] = glob[k]
    end
end

"Transfer - vector, non-trivial map"
function transfer_to_local!(loc::AbstractVector, glob::AbstractVector, M::Jutul.FiniteVolumeGlobalMap{R}) where R
    for c in eachindex(loc)
        @inbounds gc = Jutul.global_cell(c, M)::R
        @inbounds loc[c] = glob[gc]
    end
end

"Transfer - matrix, non-trivial map"
function transfer_to_local!(loc::AbstractMatrix, glob::AbstractMatrix, M::Jutul.FiniteVolumeGlobalMap{R}) where R
    nd, nc = size(loc)
    for c in 1:nc
        @inbounds gc = Jutul.global_cell(c, M)::R
        for d in 1:nd
            @inbounds loc[d, c] = glob[d, gc]
        end
    end
end
