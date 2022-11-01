struct PhaseViscosities <: PhaseVariables end
Jutul.default_value(model, v::PhaseViscosities) = 1e-3
Jutul.minimum_value(v::PhaseViscosities) = eps(Float64)
