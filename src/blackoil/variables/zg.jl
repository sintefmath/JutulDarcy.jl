function update_gas_fraction!(zg, state, zgvar, model, dx)
    abs_max_base = absolute_increment_limit(zgvar)
    maxval, minval = maximum_value(zgvar), minimum_value(zgvar)
    active_cells = active_entities(model.domain, Cells())
    sys = model.system
    F_rs = sys.saturation_table
    rhoOS = sys.rhoLS
    rhoGS = sys.rhoVS

    pressure = state.Pressure
    sat_chop = zgvar.sat_chop
    update_gas_mass_fractions_inner!(zg, dx, active_cells, pressure, sat_chop, rhoOS, rhoGS, F_rs, abs_max_base, maxval, minval)
end

function update_gas_mass_fractions_inner!(zg, dx, active_cells, pressure, sat_chop, rhoOS, rhoGS, F_rs, abs_max_base, maxval, minval)
    F = eltype(dx)
    ϵ = 1e-4
    @inbounds for (i, cell) in enumerate(active_cells)
        z = value(zg[cell])::F
        dz = dx[i]
        abs_max = abs_max_base
        if sat_chop
            rs_sat = F_rs(value(pressure[cell]))
            z_max = max_dissolved_gas_fraction(rs_sat, rhoOS, rhoGS)
            if (dz > 0 && z < z_max) || (dz < 0 && z > z_max)
                # We might cross the bubble point. Account for this as absolute limit
                abs_max_sat = abs(z_max - z) + ϵ
            else
                abs_max_sat = Inf
            end
            abs_max = min(abs_max, abs_max_sat)
        end
        dz = Jutul.choose_increment(z, dz, abs_max, nothing, minval, maxval)
        zg[cell] += dz
    end
end

@jutul_secondary function update_as_secondary!(rs, m::Rs, model::BlackOilModelGasFraction, param, PhaseState, Pressure, GasMassFraction)
    tab = model.system.saturation_table
    rhoS = param[:reference_densities]
    rhoOS = rhoS[2]
    rhoGS = rhoS[3]
    @inbounds for i in eachindex(PhaseState)
        p = Pressure[i]
        phase_state = PhaseState[i]
        if phase_state == GasOnly
            rs_i = 0
        else
            if phase_state == OilAndGas
                z_g = max_dissolved_gas_fraction(tab(p), rhoOS, rhoGS)
            else
                z_g = GasMassFraction[i]
            end
            rs_i = rhoOS*z_g/(rhoGS*(1-z_g))
        end
        rs[i] = rs_i
    end
end

max_dissolved_gas_fraction(rs, rhoOS, rhoGS) = rs*rhoGS/(rhoOS + rs*rhoGS)

@inline oil_saturation(zo, rsSat, rhoOS, rhoGS, bO, bG) = zo*rhoGS*bG/(rhoOS*bO + zo*(rhoGS*bG - rhoOS*bO - rhoGS*bO*rsSat))

@jutul_secondary function update_as_secondary!(s, SAT::Saturations, model::BlackOilModelGasFraction, param, ImmiscibleSaturation, PhaseState, GasMassFraction, ShrinkageFactors, Rs)
    # tb = minbatch(model.context)
    nph, nc = size(s)
    a, l, v = 1, 2, 3
    rhoS = param[:reference_densities]
    rhoOS = rhoS[l]
    rhoGS = rhoS[v]
    @inbounds for i = 1:nc
        sw = ImmiscibleSaturation[i]
        s[a, i] = sw
        @inbounds if PhaseState[i] == OilAndGas || PhaseState[i] == GasOnly
            rs = Rs[i]
            bO = ShrinkageFactors[l, i]
            bG = ShrinkageFactors[v, i]

            zo = 1 - GasMassFraction[i]
            so = oil_saturation(zo, rs, rhoOS, rhoGS, bO, bG)
            s_og = 1 - sw
            s[l, i] = s_og*so
            s[v, i] = s_og*(1 - so)
        else # OilOnly
            s[l, i] = 1 - sw
            s[v, i] = 0
        end
    end
end

@jutul_secondary function update_as_secondary!(phase_state, m::BlackOilPhaseState, model::BlackOilModelGasFraction, param, Pressure, GasMassFraction)
    tab = model.system.saturation_table
    rhoS = param[:reference_densities]
    rhoOS = rhoS[2]
    rhoGS = rhoS[3]
    @inbounds for i in eachindex(phase_state)
        p = Pressure[i]
        z_g = GasMassFraction[i]
        z_g_bub = max_dissolved_gas_fraction(tab(p), rhoOS, rhoGS)
        if value(z_g) == 1
            new_state = GasOnly
        elseif z_g_bub < z_g
            new_state = OilAndGas
        else
            new_state = OilOnly
        end
        phase_state[i] = new_state
    end
end

