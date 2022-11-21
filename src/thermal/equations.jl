function select_equations!(eqs, sys::ThermalSystem, model)
    disc = model.domain.discretizations.heat_flow
    eqs[:energy_conservation] = ConservationLaw(disc, :TotalThermalEnergy, 1)
end
