function reset_nldd_logger!(log, active, forces, dt)
    log.active = active
    log.dt = dt
    old_forces = log.forces
    log.forces_changed = forces != old_forces && !isnothing(old_forces)
    log.forces = forces
    empty!(log.reports)
end

function store_metadata_nldd_logger!(log, report, active)
    log.active = active
    push!(log.reports, report)
end
