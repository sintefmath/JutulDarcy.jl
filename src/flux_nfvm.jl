function flux_primitives(face, state, model, flux_type::Jutul.DefaultFlux, tpfa::Jutul.NFVM.NFVMDiscretization, upw)
    return nothing
end

function kgrad_common(face, state, model, tpfa::Jutul.NFVM.NFVMDiscretization)
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
        z = state.AdjustedCellDepths
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

@inline function JutulDarcy.gradient(X::AbstractVector, hf::Jutul.NFVM.NFVMDiscretization)
    l, r = Jutul.NFVM.cell_pair(hf)
    return @inbounds X[r] - X[l]
end

@inline function JutulDarcy.gradient(X::AbstractMatrix, i, hf::Jutul.NFVM.NFVMDiscretization)
    l, r = Jutul.NFVM.cell_pair(hf)
    return @inbounds X[i, r] - X[i, l]
end

function Jutul.get_dependencies(pot::PhasePotentials, model)
    deps = Symbol[:Pressure, :AdjustedCellDepths, :PhaseMassDensities]
    has_pc = !isnothing(get_variable(model, :CapillaryPressure, throw = false))
    if has_pc
        push!(deps, :CapillaryPressure)
    end
    return tuple(deps...)
end

function update_secondary_variable!(pot, pot_def::PhasePotentials, model, state, ix = entity_eachindex(kr))
    pc, ref_index = capillary_pressure(model, state)
    AdjustedCellDepths = state.AdjustedCellDepths
    Pressure = state.Pressure
    PhaseMassDensities = state.PhaseMassDensities
    g = gravity_constant
    for i in ix
        for ph in axes(pot, 1)
            z = AdjustedCellDepths[i]
            p = Pressure[i]
            rho = PhaseMassDensities[ph, i]
            if isnothing(pc) || ref_index == ph
                p_ph = p
            else
                pos = ph - (ph > ref_index)
                p_ph = p + pc[pos, i]
            end
            pot[ph, i] = p_ph - g*z*rho
        end
    end
    return pot
end

function Jutul.default_values(model, ::AdjustedCellDepths)
    data_domain = model.data_domain
    nc = number_of_cells(data_domain)
    v = missing
    if haskey(data_domain, :cell_centroids)
        fc = data_domain[:cell_centroids]
        if size(fc, 1) == 3
            # Make sure that this depth is always negative
            z = fc[3, :]
            zmax = maximum(z)
            v = @. z - zmax
            @assert all(x -> x <= 0, v)
        end
    end
    if ismissing(v)
        v = zeros(nc)
    end
    return v
end
