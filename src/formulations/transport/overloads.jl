function Jutul.select_primary_variables!(vars, ::TransportFormulation, model::TransportModel)
    delete!(vars, :Pressure)
    vars[:TotalSaturation] = TotalSaturation()
end

function Jutul.select_parameters!(vars, ::TransportFormulation, model::TransportModel)
    vars[:Pressure] = Pressure()
    vars[:TotalVolumetricFlux] = TotalVolumetricFlux()
end
