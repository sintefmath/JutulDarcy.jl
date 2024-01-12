# Total masses
@jutul_secondary function update_total_masses!(totmass, tv::TotalMasses, model::StandardBlackOilModel,
                                                                                                    Rs,
                                                                                                    Rv,
                                                                                                    ShrinkageFactors,
                                                                                                    PhaseMassDensities,
                                                                                                    Saturations,
                                                                                                    FluidVolume,
                                                                                                    ix)
    sys = model.system
    rhoS = reference_densities(sys)
    ind = phase_indices(sys)
    has_water = has_other_phase(sys)
    for cell in ix
        @inbounds @views blackoil_mass!(
            totmass[:, cell],
            FluidVolume, PhaseMassDensities, Rs, Rv, ShrinkageFactors, Saturations,
            rhoS, cell, ind, has_water)
    end
end

@jutul_secondary function update_total_masses!(totmass, tv::TotalMasses, model::VapoilBlackOilModel,
                                                                                                    Rv,
                                                                                                    ShrinkageFactors,
                                                                                                    PhaseMassDensities,
                                                                                                    Saturations,
                                                                                                    FluidVolume,
                                                                                                    ix)
    sys = model.system
    rhoS = reference_densities(sys)
    ind = phase_indices(sys)
    has_water = has_other_phase(sys)
    for cell in ix
        @inbounds @views blackoil_mass!(
            totmass[:, cell],
            FluidVolume, PhaseMassDensities, missing, Rv, ShrinkageFactors, Saturations,
            rhoS, cell, ind, has_water)
    end
end

@jutul_secondary function update_total_masses!(totmass, tv::TotalMasses, model::DisgasBlackOilModel,
                                                                                                    Rs,
                                                                                                    ShrinkageFactors,
                                                                                                    PhaseMassDensities,
                                                                                                    Saturations,
                                                                                                    FluidVolume,
                                                                                                    ix)
    sys = model.system
    rhoS = reference_densities(sys)
    ind = phase_indices(sys)
    has_water = has_other_phase(sys)
    for cell in ix
        @inbounds @views blackoil_mass!(
            totmass[:, cell],
            FluidVolume, PhaseMassDensities, Rs, missing, ShrinkageFactors, Saturations,
            rhoS, cell, ind, has_water)
    end
end

Base.@propagate_inbounds @inline function blackoil_mass!(M, pv, ρ, Rs, Rv, b, S, rhoS, cell, phase_indices, has_water)
    has_disgas = !ismissing(Rs)
    has_vapoil = !ismissing(Rv)
    Φ = pv[cell]
    if has_water
        a, l, v = phase_indices
        # Water is trivial
        M[a] = Φ*ρ[a, cell]*S[a, cell]
    else
        l, v = phase_indices
    end
    bO = b[l, cell]
    bG = b[v, cell]
    sO = S[l, cell]
    sG = S[v, cell]

    if has_vapoil
        # Oil is in both phases
        rv = Rv[cell]
        M_o = bO*sO + bG*sG*rv
    else
        M_o = bO*sO
    end
    M[l] = Φ*rhoS[l]*M_o
    if has_disgas
        # Gas is in both phases
        rs = Rs[cell]
        M_g = bG*sG + bO*sO*rs
    else
        M_g = bG*sG
    end
    M[v] = Φ*rhoS[v]*M_g
end
