struct TotalSaturation <: Jutul.ScalarVariable end
Jutul.minimum_value(::TotalSaturation) = 1e-10

struct TotalVolumetricFlux <: Jutul.ScalarVariable end
Jutul.associated_entity(::TotalVolumetricFlux) = Faces()

struct TotalSaturationCorrectedVariable <: PhaseVariables end

@jutul_secondary function total_sat_correction(v, ts::TotalSaturationCorrectedVariable, model::TransportModel, UncorrectedVariable, TotalSaturation, ix)
    for i in ix
        sT = TotalSaturation[i]
        for j in axes(v, 1)
            v[j, i] = sT*UncorrectedVariable[j, i]
        end
    end
end