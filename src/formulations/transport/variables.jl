struct TotalSaturation <: Jutul.ScalarVariable end

Jutul.minimum_value(::TotalSaturation) = 1e-10
