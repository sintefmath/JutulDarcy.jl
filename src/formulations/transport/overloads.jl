function Jutul.select_primary_variables!(vars, ::TransportFormulation, model::TransportModel)
    delete!(vars, :Pressure)
    vars[:TotalSaturation] = TotalSaturation()
end

function Jutul.select_parameters!(vars, ::TransportFormulation, model::TransportModel)
    vars[:Pressure] = Pressure()
    vars[:TotalVolumetricFlux] = TotalVolumetricFlux()
end

function Jutul.select_equations!(eqs, ::TransportFormulation, model::TransportModel)
    lbl = :mass_conservation
    @assert haskey(eqs, lbl)
    eq = eqs[lbl]

    s = conserved_symbol(eq)
    @assert s == :TotalMasses
    N = number_of_equations_per_entity(model, eq)
    eqs[lbl] = ConservationLaw(eq.flow_discretization, s, N, flux = TotalSaturationFlux())
end