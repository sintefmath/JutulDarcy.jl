function heat_flux_boundary_condition(cells, domain, heat_flux)
    bc = []
    heat_flux_boundary_condition!(bc, domain, cells, heat_flux)
    return [i for i in bc]
end

function heat_flux_boundary_condition!(bc, domain, cells, heat_flux)
    n = length(cells)
    if heat_flux isa Real
        heat_flux = fill(heat_flux, n)
    end
    for (cell, q) in zip(cells, heat_flux)
        bc_c = HeatFluxBoundaryCondition(domain, cell, q)
        push!(bc, bc_c)
    end

    return bc
end

function apply_force_to_equation!(acc, storage, model::SimulationModel{D, S}, eq::ConservationLaw{:TotalMasses}, eq_s, bc::V, time) where {V <: HeatFluxBoundaryCondition, D, S<:MultiPhaseSystem}

end

function apply_force_to_equation!(acc, storage, model::SimulationModel{D, S}, eq::ConservationLaw{:TotalThermalEnergy}, eq_s, bc::V, time) where {V <: HeatFluxBoundaryCondition, D, S<:MultiPhaseSystem}
    c = bc.cell
    acc_i = view(acc, :, c)
    acc_i .-= bc.heat_flux
end