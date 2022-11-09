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

struct BlackOilPhaseState <: ScalarVariable

end

Jutul.default_value(model, ::BlackOilPhaseState) = OilAndGas
Jutul.initialize_secondary_variable_ad!(state, model, var::BlackOilPhaseState, arg...; kwarg...) = state

struct Rs <: ScalarVariable end
struct Rv <: ScalarVariable end

Base.@kwdef struct BlackOilUnknown{R} <: ScalarVariable
    dr_max::R = Inf
    ds_max::R = 0.2
end

export BlackOilX
struct BlackOilX{T}
    val::T
    phases_present::PresentPhasesBlackOil
    sat_close::Bool
    function BlackOilX(val::T, phases = OilAndGas, sat_close = false) where T<:Real
        return new{T}(val, phases, sat_close)
    end
end

Jutul.default_value(model, ::BlackOilUnknown) = (NaN, OilAndGas, false) # NaN, Oil+Gas, away from bubble point
function Jutul.initialize_primary_variable_ad!(state, model, pvar::BlackOilUnknown, symb, npartials; offset, kwarg...)
    pre = state[symb]
    vals = map(x -> x.val, pre)
    ad_vals = allocate_array_ad(vals, diag_pos = offset + 1, context = model.context, npartials = npartials; kwarg...)
    state[symb] = map((v, x) -> BlackOilX(v, x.phases_present, false), ad_vals, pre)
    return state
end

@jutul_secondary function update_as_secondary!(b, ρ::DeckShrinkageFactors, model::DisgasBlackOilModel, Pressure, Rs)
    pvt, reg = ρ.pvt, ρ.regions
    nph, nc = size(b)
    tb = minbatch(model.context, nc)

    w, o, g = phase_indices(model.system)
    bO = pvt[o]
    bG = pvt[g]
    bW = pvt[w]
    @inbounds @batch minbatch = tb for i in 1:nc
        p = Pressure[i]
        rs = Rs[i]
        b[w, i] = shrinkage(bW, reg, p, i)
        b[o, i] = shrinkage(bO, reg, p, rs, i)
        b[g, i] = shrinkage(bG, reg, p, i)
    end
end

# Shrinkage factors for all three cases
@jutul_secondary function update_as_secondary!(b, ρ::DeckShrinkageFactors, model::StandardBlackOilModel, Pressure, Rs, Rv)
    pvt, reg = ρ.pvt, ρ.regions
    nph, nc = size(b)
    tb = minbatch(model.context, nc)

    w, o, g = phase_indices(model.system)
    bO = pvt[o]
    bG = pvt[g]
    bW = pvt[w]
    @inbounds @batch minbatch = tb for i in 1:nc
        p = Pressure[i]
        rv = Rv[i]
        rs = Rs[i]
        b[w, i] = shrinkage(bW, reg, p, i)
        b[o, i] = shrinkage(bO, reg, p, rs, i)
        b[g, i] = shrinkage(bG, reg, p, rv, i)
    end
end

@jutul_secondary function update_as_secondary!(b, ρ::DeckShrinkageFactors, model::VapoilBlackOilModel, Pressure, Rv)
    pvt, reg = ρ.pvt, ρ.regions
    nph, nc = size(b)
    tb = minbatch(model.context, nc)

    w, o, g = phase_indices(model.system)
    bO = pvt[o]
    bG = pvt[g]
    bW = pvt[w]
    @inbounds @batch minbatch = tb for i in 1:nc
        p = Pressure[i]
        rv = Rv[i]
        b[w, i] = shrinkage(bW, reg, p, i)
        b[o, i] = shrinkage(bO, reg, p, i)
        b[g, i] = shrinkage(bG, reg, p, rv, i)
    end
end

@jutul_secondary function update_as_secondary!(μ, ρ::DeckViscosity, model::DisgasBlackOilModel, Pressure, Rs)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    nph, nc = size(μ)

    w, o, g = phase_indices(model.system)
    muW = pvt[w]
    muO = pvt[o]
    muG = pvt[g]

    mb = minbatch(model.context, nc)
    @inbounds @batch minbatch = mb for i = 1:nc
        p = Pressure[i]
        rs = Rs[i]
        μ[w, i] = viscosity(muW, reg, p, i)
        μ[o, i] = viscosity(muO, reg, p, rs, i)
        μ[g, i] = viscosity(muG, reg, p, i)
    end
end

# Density for all three cases
@jutul_secondary function update_as_secondary!(rho, m::DeckDensity, model::StandardBlackOilModel, Rs, Rv, ShrinkageFactors)
    b = ShrinkageFactors
    sys = model.system
    w, o, g = phase_indices(sys)
    rhoWS, rhoOS, rhoGS = reference_densities(sys)
    n = size(rho, 2)
    mb = minbatch(model.context, n)
    @inbounds @batch minbatch = mb for i = 1:n
        rho[w, i] = b[w, i]*rhoWS
        rho[o, i] = b[o, i]*(rhoOS + Rs[i]*rhoGS)
        rho[g, i] = b[g, i]*(rhoGS + Rv[i]*rhoOS)
    end
end

@jutul_secondary function update_as_secondary!(rho, m::DeckDensity, model::VapoilBlackOilModel, Rv, ShrinkageFactors)
    b = ShrinkageFactors
    sys = model.system
    w, o, g = phase_indices(sys)
    rhoWS, rhoOS, rhoGS = reference_densities(sys)
    n = size(rho, 2)
    mb = minbatch(model.context, n)
    @inbounds @batch minbatch = mb for i = 1:n
        rho[w, i] = b[w, i]*rhoWS
        rho[o, i] = b[o, i]*rhoOS
        rho[g, i] = b[g, i]*(rhoGS + Rv[i]*rhoOS)
    end
end

@jutul_secondary function update_as_secondary!(rho, m::DeckDensity, model::DisgasBlackOilModel, Rs, ShrinkageFactors)
    b = ShrinkageFactors
    sys = model.system
    w, o, g = phase_indices(sys)
    rhoWS, rhoOS, rhoGS = reference_densities(sys)
    n = size(rho, 2)
    mb = minbatch(model.context, n)
    @inbounds @batch minbatch = mb for i = 1:n
        rho[w, i] = b[w, i]*rhoWS
        rho[o, i] = b[o, i]*(rhoOS + Rs[i]*rhoGS)
        rho[g, i] = b[g, i]*rhoGS
    end
end

# Total masses
@jutul_secondary function update_as_secondary!(totmass, tv::TotalMasses, model::StandardBlackOilModel,
                                                                                                    Rs,
                                                                                                    Rv,
                                                                                                    ShrinkageFactors,
                                                                                                    PhaseMassDensities,
                                                                                                    Saturations,
                                                                                                    FluidVolume)
    sys = model.system
    rhoS = reference_densities(sys)
    ind = phase_indices(sys)
    nc = size(totmass, 2)
    tb = minbatch(model.context, nc)
    @batch minbatch = tb for cell = 1:nc
        @inbounds @views blackoil_mass!(totmass[:, cell], FluidVolume, PhaseMassDensities, Rs, Rv, ShrinkageFactors, Saturations, rhoS, cell, ind)
    end
end

Base.@propagate_inbounds function blackoil_mass!(M, pv, ρ, Rs, Rv, b, S, rhoS, cell, phase_indices)
    a, l, v = phase_indices
    bO = b[l, cell]
    bG = b[v, cell]
    rs = Rs[cell]
    rv = Rv[cell]
    sO = S[l, cell]
    sG = S[v, cell]
    Φ = pv[cell]

    # Water is trivial
    M[a] = Φ*ρ[a, cell]*S[a, cell]
    # Oil is only in oil phase
    M[l] = Φ*rhoS[l]*(bO*sO + bG*sG*rv)
    # Gas is in both phases
    M[v] = Φ*rhoS[v]*(bG*sG + bO*sO*rs)
end

@jutul_secondary function update_as_secondary!(totmass, tv::TotalMasses, model::StandardBlackOilModel,
                                                                                                    Rv,
                                                                                                    ShrinkageFactors,
                                                                                                    PhaseMassDensities,
                                                                                                    Saturations,
                                                                                                    FluidVolume)
    sys = model.system
    rhoS = reference_densities(sys)
    ind = phase_indices(sys)
    nc = size(totmass, 2)
    tb = minbatch(model.context, nc)
    @batch minbatch = tb for cell = 1:nc
        @inbounds @views blackoil_vapoil!(totmass[:, cell], FluidVolume, PhaseMassDensities, Rv, ShrinkageFactors, Saturations, rhoS, cell, ind)
    end
end

Base.@propagate_inbounds function blackoil_vapoil!(M, pv, ρ, Rv, b, S, rhoS, cell, phase_indices)
    a, l, v = phase_indices
    bO = b[l, cell]
    bG = b[v, cell]
    rv = Rv[cell]
    sO = S[l, cell]
    sG = S[v, cell]
    Φ = pv[cell]

    # Water is trivial
    M[a] = Φ*ρ[a, cell]*S[a, cell]
    # Oil is only in oil phase
    M[l] = Φ*rhoS[l]*(bO*sO + bG*sG*rv)
    # Gas is in both phases
    M[v] = Φ*rhoS[v]*bG*sG
end

@jutul_secondary function update_as_secondary!(totmass, tv::TotalMasses, model::DisgasBlackOilModel,
                                                                                                    Rs,
                                                                                                    ShrinkageFactors,
                                                                                                    PhaseMassDensities,
                                                                                                    Saturations,
                                                                                                    FluidVolume)
    sys = model.system
    rhoS = reference_densities(sys)
    ind = phase_indices(sys)
    nc = size(totmass, 2)
    tb = minbatch(model.context, nc)
    @batch minbatch = tb for cell = 1:nc
        @inbounds @views blackoil_mass_disgas!(totmass[:, cell], FluidVolume, PhaseMassDensities, Rs, ShrinkageFactors, Saturations, rhoS, cell, ind)
    end
end

Base.@propagate_inbounds function blackoil_mass_disgas!(M, pv, ρ, Rs, b, S, rhoS, cell, phase_indices)
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

struct SurfaceVolumeMobilities <: PhaseVariables end

@jutul_secondary function update_as_secondary!(b_mob, var::SurfaceVolumeMobilities, model,
                                                        ShrinkageFactors,
                                                        PhaseViscosities,
                                                        RelativePermeabilities)
    # For blackoil, the main upwind term
    mb = minbatch(model.context)
    @batch minbatch = mb for i in axes(b_mob, 2)
        @inbounds for ph in axes(b_mob, 1)
            b_mob[ph, i] = ShrinkageFactors[ph, i]*RelativePermeabilities[ph, i]/PhaseViscosities[ph, i]
        end
    end
end

include("zg.jl")
include("varswitch.jl")
