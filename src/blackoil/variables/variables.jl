Base.@kwdef struct GasMassFraction{R} <: ScalarVariable
    dz_max::R = 0.2
    sat_chop::Bool = true
end

maximum_value(::GasMassFraction) = 1.0
minimum_value(::GasMassFraction) = 1e-12
absolute_increment_limit(s::GasMassFraction) = s.dz_max

function update_primary_variable!(state, p::GasMassFraction, state_symbol, model, dx)
    s = state[state_symbol]
    update_gas_fraction!(s, state, p, model, dx)
end

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

struct BlackOilPhaseState <: ScalarVariable

end

Jutul.default_value(model, ::BlackOilPhaseState) = OilAndGas
Jutul.initialize_secondary_variable_ad!(state, model, var::BlackOilPhaseState, arg...; kwarg...) = state

max_dissolved_gas_fraction(rs, rhoOS, rhoGS) = rs*rhoGS/(rhoOS + rs*rhoGS)

@jutul_secondary function update_as_secondary!(phase_state, m::BlackOilPhaseState, model::SimulationModel{D, S}, param, Pressure, GasMassFraction) where {D, S<:BlackOilSystem}
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

struct Rs <: ScalarVariable end

@jutul_secondary function update_as_secondary!(rs, m::Rs, model::SimulationModel{D, S}, param, PhaseState, Pressure, GasMassFraction) where {D, S<:BlackOilSystem}
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

@jutul_secondary function update_as_secondary!(b, ρ::DeckShrinkageFactors, model::SimulationModel{D, StandardBlackOilSystem{T, true, R}}, param, Pressure, Rs) where {D, T, R}
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    tb = minbatch(model.context)
    nph, nc = size(b)

    w = 1
    g = 3
    o = 2
    bO = pvt[o]
    bG = pvt[g]
    bW = pvt[w]
    @batch minbatch = tb for i in 1:nc
        @inbounds p = Pressure[i]
        @inbounds rs = Rs[i]
        b[w, i] = shrinkage(bW, reg, p, i)
        b[o, i] = shrinkage(bO, reg, p, rs, i)
        b[g, i] = shrinkage(bG, reg, p, i)
    end
end

@jutul_secondary function update_as_secondary!(μ, ρ::DeckViscosity, model::SimulationModel{D, StandardBlackOilSystem{T, true, R}}, param, Pressure, Rs) where {D, T, R}
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    tb = minbatch(model.context)
    nph, nc = size(μ)

    w = 1
    o = 2
    g = 3
    muW = pvt[w]
    muO = pvt[o]
    muG = pvt[g]

    # @batch minbatch = tb for i in 1:nc
    @inbounds for i = 1:nc
        p = Pressure[i]
        rs = Rs[i]
        μ[w, i] = viscosity(muW, reg, p, i)
        μ[o, i] = viscosity(muO, reg, p, rs, i)
        μ[g, i] = viscosity(muG, reg, p, i)
    end
end

@jutul_secondary function update_as_secondary!(rho, m::DeckDensity, model::SimulationModel{D, S}, param, Rs, ShrinkageFactors) where {D, S<:BlackOilSystem}
    # sys = model.system
    b = ShrinkageFactors
    rhoS = param[:reference_densities]
    rhoWS = rhoS[1]
    rhoOS = rhoS[2]
    rhoGS = rhoS[3]
    # pvt, reg = ρ.pvt, ρ.regions
    # eos = sys.equation_of_state
    w = 1
    o = 2
    g = 3
    n = size(rho, 2)
    @inbounds for i = 1:n
        rho[w, i] = b[w, i]*rhoWS
        rho[o, i] = b[o, i]*(rhoOS + Rs[i]*rhoGS)
        rho[g, i] = b[g, i]*rhoGS
    end
end

@jutul_secondary function update_as_secondary!(totmass, tv::TotalMasses, model::SimulationModel{G,S}, param,
                                                                                                    Rs,
                                                                                                    ShrinkageFactors,
                                                                                                    PhaseMassDensities,
                                                                                                    Saturations,
                                                                                                    FluidVolume) where {G,S<:BlackOilSystem}
    rhoS = tuple(param[:reference_densities]...)
    tb = minbatch(model.context)
    sys = model.system
    nc = size(totmass, 2)
    # @batch minbatch = tb for cell = 1:nc
    for cell = 1:nc
        @inbounds @views blackoil_mass!(totmass[:, cell], FluidVolume, PhaseMassDensities, Rs, ShrinkageFactors, Saturations, rhoS, cell, (1,2,3))
    end
end

Base.@propagate_inbounds function blackoil_mass!(M, pv, ρ, Rs, b, S, rhoS, cell, phase_indices)
    a, l, v = phase_indices
    bO = b[l, cell]
    bG = b[v, cell]
    rs = Rs[cell]
    sO = S[l, cell]
    sG = S[v, cell]
    Φ = pv[cell]

    # Water is trivial
    M[a] = Φ*ρ[a, cell]*S[a, cell]
    # Oil is only in oil phase
    M[l] = Φ*rhoS[l]*bO*sO
    # Gas is in both phases
    M[v] = Φ*rhoS[v]*(bG*sG + bO*sO*rs)
end

@inline oil_saturation(zo, rsSat, rhoOS, rhoGS, bO, bG) = zo*rhoGS*bG/(rhoOS*bO + zo*(rhoGS*bG - rhoOS*bO - rhoGS*bO*rsSat))

@jutul_secondary function update_as_secondary!(s, SAT::Saturations, model::SimulationModel{D, S}, param, ImmiscibleSaturation, PhaseState, GasMassFraction, ShrinkageFactors, Rs) where {D, S<:BlackOilSystem}
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
