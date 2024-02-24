@jutul_secondary function update_fluid_internal_energy!(U, fe::FluidInternalEnergy, model::ThermalImmiscibleModel, Temperature, ComponentHeatCapacity, ix)
    @assert size(U) == size(ComponentHeatCapacity) "This fluid internal energy implementation assumes immiscible phases."
    for i in ix
        T = Temperature[i]
        for c in axes(U, 1)
            U[c, i] = ComponentHeatCapacity[c, i]*T
        end
    end
end

@jutul_secondary function update_fluid_enthalpy!(H, fe::FluidEnthalpy, model::ThermalModel, FluidInternalEnergy, Pressure, PhaseMassDensities, ix)
    for i in ix
        p = Pressure[i]
        for ph in axes(H, 1)
            H[ph, i] = FluidInternalEnergy[ph, i] + p/PhaseMassDensities[ph, i]
        end
    end
end

@jutul_secondary function update_rock_internal_energy!(U_r, e::RockInternalEnergy, model::ThermalModel, RockHeatCapacity, Temperature, ix)
    for i in ix
        U_r[i] = RockHeatCapacity[i]*Temperature[i]
    end
end

@jutul_secondary function update_total_thermal_energy!(E_total, te::TotalThermalEnergy, model::ThermalModel, Saturations, PhaseMassDensities, FluidInternalEnergy, RockDensity, RockInternalEnergy, BulkVolume, FluidVolume, ix)
    U_f = FluidInternalEnergy
    U_r = RockInternalEnergy
    ρ_f = PhaseMassDensities
    S = Saturations
    V = BulkVolume
    for i in ix
        V_f = FluidVolume[i]
        E_i = RockDensity[i]*U_r[i]*(V[i] - V_f)
        for ph in axes(U_f, 1)
            E_i += ρ_f[ph, i]*S[ph, i]*U_f[ph, i]*V_f
        end
        E_total[i] = E_i
    end
end
