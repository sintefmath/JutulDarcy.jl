@jutul_secondary function update_fluid_internal_energy!(U, fe::FluidInternalEnergy, model, Temperature, ComponentHeatCapacity, ix)
    C = ComponentHeatCapacity
    @assert size(U) == size(C) "This fluid internal energy implementation assumes immiscible phases."
    for i in ix
        T = Temperature[i]
        for c in axes(U, 1)
            U[c, i] = C[c, i]*T
        end
    end
end

@jutul_secondary function update_fluid_internal_energy!(U, fe::FluidInternalEnergy, model::SinglePhaseModel, Temperature, ComponentHeatCapacity, ix)
    C = ComponentHeatCapacity
    @assert size(U) == size(C) "This fluid internal energy implementation assumes immiscible phases."
    for i in ix
        T = Temperature[i]
        for c in axes(U, 1)
            U[c, i] = C[c, i]*T
        end
    end
end

@jutul_secondary function update_fluid_internal_energy!(U, fe::FluidInternalEnergy, model::CompositionalModel, Temperature, ComponentHeatCapacity, LiquidMassFractions, VaporMassFractions, ix)
    fsys = flow_system(model.system)
    C = ComponentHeatCapacity
    X = LiquidMassFractions
    Y = VaporMassFractions
    if has_other_phase(fsys)
        a, l, v = phase_indices(fsys)
        offset = 1
    else
        l, v = phase_indices(fsys)
        offset = 0
    end
    for i in ix
        T = Temperature[i]
        if has_other_phase(fsys)
            C_w = C[1, i]
            U[a, i] = C_w*T
        end
        U_v = 0.0
        U_l = 0.0
        for c in axes(X, 1)
            C_i = C[c+offset, i]
            U_l += C_i*X[c, i]
            U_v += C_i*Y[c, i]
        end
        U[l, i] = U_l*T
        U[v, i] = U_v*T
    end
end

@jutul_secondary function update_fluid_enthalpy!(H, fe::FluidEnthalpy, model, FluidInternalEnergy, Pressure, PhaseMassDensities, ix)
    for i in ix
        p = Pressure[i]
        for ph in axes(H, 1)
            H[ph, i] = FluidInternalEnergy[ph, i] + p/PhaseMassDensities[ph, i]
        end
    end
end

@jutul_secondary function update_rock_internal_energy!(U_r, e::RockInternalEnergy, model, RockHeatCapacity, Temperature, ix)
    for i in ix
        U_r[i] = RockHeatCapacity[i]*Temperature[i]
    end
end

@jutul_secondary function update_total_thermal_energy!(E_total, te::TotalThermalEnergy, model, Saturations, PhaseMassDensities, FluidInternalEnergy, RockDensity, RockInternalEnergy, BulkVolume, FluidVolume, ix)
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

const MSWellDomain = DiscretizedDomain{<:MultiSegmentWell}
const MSWellFlowModel = SimulationModel{<:MSWellDomain, <:MultiPhaseSystem}
@jutul_secondary function update_total_thermal_energy!(E_total, te::TotalThermalEnergy, model::MSWellFlowModel, Saturations, PhaseMassDensities, FluidInternalEnergy, MaterialDensities, MaterialInternalEnergy, BulkVolume, FluidVolume, ix)
    update_total_thermal_energy!(E_total, te::TotalThermalEnergy, nothing, Saturations, PhaseMassDensities, FluidInternalEnergy, MaterialDensities, MaterialInternalEnergy, BulkVolume, FluidVolume, ix)
end

@jutul_secondary function update_material_internal_energy!(U_m, e::MaterialInternalEnergy, model::MSWellFlowModel, MaterialHeatCapacities, Temperature, ix)
    for i in ix
        U_m[i] = MaterialHeatCapacities[i]*Temperature[i]
    end
end