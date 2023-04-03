struct TotalSaturation <: Jutul.ScalarVariable end
Jutul.minimum_value(::TotalSaturation) = 1e-10

struct TotalVolumetricFlux <: Jutul.ScalarVariable end
Jutul.associated_entity(::TotalVolumetricFlux) = Faces()
