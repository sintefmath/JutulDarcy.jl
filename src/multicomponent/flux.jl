@inline function component_mass_fluxes!(q, face, state, model::SimulationModel{<:Any, <:CompositionalSystem, <:Any, <:Any}, kgrad, upw)
    sys = model.system
    aqua = Val(has_other_phase(sys))
    ph_ix = phase_indices(sys)

    X = state.LiquidMassFractions
    Y = state.VaporMassFractions
    kdisc = kgrad_common(face, state, model, kgrad)
    q = compositional_fluxes!(q, face, state,  X, Y, model, kgrad, kdisc, upw, aqua, ph_ix)
    return q
end

@inline function compositional_fluxes!(q, face, state, X, Y, model, kgrad, kdisc, upw, aqua::Val{false}, phase_ix)
    nc = size(X, 1)
    l, v = phase_ix
    q_l = darcy_phase_mass_flux(face, l, state, model, kgrad, upw, kdisc)
    q_v = darcy_phase_mass_flux(face, v, state, model, kgrad, upw, kdisc)

    q = inner_compositional!(q, X, Y, q_l, q_v, upw, nc)
    return q
end

@inline function compositional_fluxes!(q, face, state, X, Y, model, kgrad, kdisc, upw, aqua::Val{true}, phase_ix)
    nc = size(X, 1)
    a, l, v = phase_ix
    q_a = darcy_phase_mass_flux(face, a, state, model, kgrad, upw, kdisc)
    q_l = darcy_phase_mass_flux(face, l, state, model, kgrad, upw, kdisc)
    q_v = darcy_phase_mass_flux(face, v, state, model, kgrad, upw, kdisc)

    q = inner_compositional!(q, X, Y, q_l, q_v, upw, nc)
    q = setindex(q, q_a, nc+1)
    return q
end

@inline function inner_compositional!(q, X, Y, q_l, q_v, upw, nc)
    for i in 1:nc
        X_f = upwind(upw, cell -> @inbounds(X[i, cell]), q_l)
        Y_f = upwind(upw, cell -> @inbounds(Y[i, cell]), q_v)

        q_i = q_l*X_f + q_v*Y_f
        q = setindex(q, q_i, i)
    end
    return q
end
