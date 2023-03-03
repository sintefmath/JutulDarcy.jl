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
    for cell in ix
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
    # Oil is in both phases
    M[l] = Φ*rhoS[l]*(bO*sO + bG*sG*rv)
    # Gas is in both phases
    M[v] = Φ*rhoS[v]*(bG*sG + bO*sO*rs)
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
    for cell in ix
        @inbounds @views blackoil_mass_vapoil!(totmass[:, cell], FluidVolume, PhaseMassDensities, Rv, ShrinkageFactors, Saturations, rhoS, cell, ind)
    end
end

Base.@propagate_inbounds function blackoil_mass_vapoil!(M, pv, ρ, Rv, b, S, rhoS, cell, phase_indices)
    a, l, v = phase_indices
    bO = b[l, cell]
    bG = b[v, cell]
    rv = Rv[cell]
    sO = S[l, cell]
    sG = S[v, cell]
    Φ = pv[cell]

    # Water is trivial
    M[a] = Φ*ρ[a, cell]*S[a, cell]
    # Oil is in both phases
    M[l] = Φ*rhoS[l]*(bO*sO + bG*sG*rv)
    # Gas is only in gas phase
    M[v] = Φ*rhoS[v]*bG*sG
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
    if has_other_phase(sys)
        for cell in ix
            @inbounds @views blackoil_mass_disgas!(totmass[:, cell], FluidVolume, PhaseMassDensities, Rs, ShrinkageFactors, Saturations, rhoS, cell, ind)
        end
    else
        for cell in ix
            @inbounds @views blackoil_mass_disgas_no_water!(totmass[:, cell], FluidVolume, PhaseMassDensities, Rs, ShrinkageFactors, Saturations, rhoS, cell, ind)
        end
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

Base.@propagate_inbounds function blackoil_mass_disgas_no_water!(M, pv, ρ, Rs, b, S, rhoS, cell, phase_indices)
    l, v = phase_indices
    bO = b[l, cell]
    bG = b[v, cell]
    rs = Rs[cell]
    sO = S[l, cell]
    sG = S[v, cell]
    Φ = pv[cell]

    # Oil is only in oil phase
    M[l] = Φ*rhoS[l]*bO*sO
    # Gas is in both phases
    M[v] = Φ*rhoS[v]*(bG*sG + bO*sO*rs)
end
