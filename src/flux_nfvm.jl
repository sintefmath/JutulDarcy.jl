function flux_primitives(face, state, model, flux_type::Jutul.DefaultFlux, tpfa::Jutul.NFVM.NFVMDiscretization, upw)
    return nothing
end

function darcy_phase_kgrad_potential(face, phase, state, model, flux_type, mpfa::D, upw, common = nothing) where {D<:Jutul.NFVM.NFVMDiscretization}
    # gΔz = tpfa.face_sign*grav[face]
    grav = state.TwoPointGravityDifference[face]
    pc, ref_index = capillary_pressure(model, state)
    p = state.Pressure
    if haskey(state, :PhasePotentials)
        g = gravity_constant
        pot = state.PhasePotentials
        dens = state.PhaseMassDensities
        z = state.CellDepths
        # grad(p - rho g z) = grad(p) - grad(rho) g z - rho g grad(z)
        # -> grad(p) - rho g grad(z) = grad(p - rho g z) + grad(rho) g z
        ∇pot = Jutul.NFVM.evaluate_flux(pot, mpfa, phase)
        ∇rho = Jutul.NFVM.evaluate_flux(dens, mpfa, phase)

        l, r = Jutul.NFVM.cell_pair(mpfa)
        z_avg = (z[l] + z[r])/2.0
        q = -(∇pot + ∇rho*g*z_avg)
    else
        # If missing potential, just do everything here.
        K∇p = Jutul.NFVM.evaluate_flux(p, mpfa, phase)
        q = -K∇p
    end
    if haskey(state, :PermeabilityMultiplier)
        K_mul = state.PermeabilityMultiplier
        m = face_average(c -> K_mul[c], tpfa)
        q *= m
    end
    return q
end

@jutul_secondary function update_phase_potentials!(pot, pot_def::PhasePotentials{false}, model, Pressure, CellDepths, PhaseMassDensities, ix)
    g = gravity_constant
    for i in ix
        for ph in axes(pot, 1)
            z = CellDepths[i]
            p = Pressure[i]
            rho = PhaseMassDensities[ph, i]
            pot[ph, i] = p - g*z*rho
        end
    end
    return pot
end

@jutul_secondary function update_phase_potentials!(pot, pot_def::PhasePotentials{true}, model, Pressure, CellDepths, PhaseMassDensities, CapillaryPressure, ix)
    error("Not supported yet")
end

function Jutul.default_values(model, ::CellDepths)
    data_domain = model.data_domain
    nc = number_of_cells(data_domain)
    v = missing
    if haskey(data_domain, :cell_centroids)
        fc = data_domain[:cell_centroids]
        if size(fc, 1) == 3
            v = fc[3, :]
        end
    end
    if ismissing(v)
        v = zeros(nc)
    end
    return v
end
