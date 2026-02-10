struct HeatFluxBoundaryCondition{I, F} <: JutulForce
    cell::I
    heat_flux::F
end

function HeatFluxBoundaryCondition(
        domain::DataDomain,
        cell::Int,
        heat_flux::Number,
    )
    return HeatFluxBoundaryCondition(cell, heat_flux)
end