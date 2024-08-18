function pressure_gradient(p, mpfa::Jutul.NFVM.NFVMDiscretization)
    return Jutul.NFVM.evaluate_flux(p, mpfa)
end

function flux_primitives(face, state, model, flux_type::Jutul.DefaultFlux, tpfa::Jutul.NFVM.NFVMDiscretization, upw)
    return nothing
end

function darcy_phase_kgrad_potential(face, phase, state, model, flux_type, mpfa::Jutul.NFVM.NFVMDiscretization, upw, common = nothing)
    # gΔz = tpfa.face_sign*grav[face]
    grav = state.TwoPointGravityDifference[face]
    pc, ref_index = capillary_pressure(model, state)
    if grav == 0.0 && isnothing(pc)
        p = state.Pressure
        K∇p = pressure_gradient(p, mpfa)
        q = -K∇p
    else
        # TODO: Sign, potential split magic, etc.
        error("Not implemented - fixme.")
        ρ_avg = face_average_density(model, state, tpfa, phase)
        Δpc = capillary_gradient(pc, l, r, phase, ref_index)
    end
    # ∇p, T_f, gΔz = common
    # l = tpfa.left
    # r = tpfa.right

    if haskey(state, :PermeabilityMultiplier)
        K_mul = state.PermeabilityMultiplier
        m = face_average(c -> K_mul[c], tpfa)
        q *= m
    end
    return q
end
