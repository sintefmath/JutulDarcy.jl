@jutul_secondary function update_pressure_factors(w, wdef::PressureReductionFactors, model::ImmiscibleModel, PhaseMassDensities, ix)
    for i in ix
        for ph in axes(w, 1)
            w[ph, i] = 1.0/PhaseMassDensities[ph, i]
        end
    end
end

function Jutul.get_dependencies(var::PressureReductionFactors, model::StandardBlackOilModel)
    sys = model.system
    deps = [:ShrinkageFactors, :FluidVolume]
    if has_disgas(sys)
        push!(deps, :Rs)
    end
    if has_vapoil(sys)
        push!(deps, :Rv)
    end
    return deps
end

function Jutul.update_secondary_variable!(w, var::PressureReductionFactors, model::StandardBlackOilModel, state, ix = entity_eachindex(V))
    sys = model.system
    T = eltype(state.Pressure)
    if has_disgas(sys)
        Rs = state.Rs
    else
        Rs = nothing
    end
    if has_vapoil(model.system)
        Rv = state.Rv
    else
        Rv = nothing
    end
    cellval(::Nothing, c) = 0.0
    cellval(x, c) = x[c]
    has_water = JutulDarcy.has_other_phase(sys)
    phase_ix = JutulDarcy.phase_indices(sys)
    rhoS = JutulDarcy.reference_densities(sys)
    rhoOS = rhoS[phase_ix.l]
    rhoGS = rhoS[phase_ix.v]
    rho = state.PhaseMassDensities
    bfactors = state.ShrinkageFactors
    for cell in ix
        pv = state.FluidVolume[cell]
        if has_water
            denw = rho[phase_ix.a, cell]
            w[phase_ix.a, cell] = 1.0/denw
        end
        bO = bfactors[phase_ix.l, cell]
        bG = bfactors[phase_ix.v, cell]
        w_o, w_g = blackoil_true_impes_og!(cellval(Rs, cell), cellval(Rv, cell), bO, bG, rhoOS, rhoGS)
        w[phase_ix.l, cell] = w_o
        w[phase_ix.v, cell] = w_g
    end
    return w
end

function blackoil_true_impes_og!(rs, rv, bO, bG, rhoOS, rhoGS,)
    α = 1.0/(1.0 - rs*rv)
    # has_free_oil = phase_state == JutulDarcy.OilAndGas || phase_state == JutulDarcy.OilOnly
    # has_free_gas = phase_state == JutulDarcy.OilAndGas || phase_state == JutulDarcy.GasOnly
    w_o = (α/rhoOS)*(1.0/bO - rs/bG)
    w_g = (α/rhoGS)*(1.0/bG - rv/bO)
    return (w_o, w_g)
end
