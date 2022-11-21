@jutul_secondary function update_fluid_internal_energy!(U, ::FluidInternalEnergy, model::ThermalModel, Saturations, PhaseMassDensities, FluidHeatCapacity, ix)
    for i in ix
        for ph in axes(U, 1)
            U[ph, i] = Saturations[ph, i]*PhaseMassDensities[ph, i]*FluidHeatCapacity[ph, i]
        end
    end
end

@jutul_secondary function update_rock_internal_energy!(U_r, ::RockInternalEnergy, model::ThermalModel, RockDensity, RockHeatCapacity, Temperature, ix)
    for i in ix
        U_r[i] = RockDensity[i]*RockHeatCapacity[i]*Temperature[i]
    end
end

@jutul_secondary function update_total_thermal_energy!(E_total, ::TotalThermalEnergy, FluidInternalEnergy, RockInternalEnergy, BulkVolume, FluidVolume)
    U_f = FluidInternalEnergy
    U_r = RockInternalEnergy
    V = BulkVolume
    for i in ix
        V_f = FluidVolume[i]
        E_i = U_r*(V[i] - V_f)
        for ph in axes(U_f, 1)
            E_i += U_f[ph, i]*V_f
        end
        E_total[i] = E_i
    end
end
