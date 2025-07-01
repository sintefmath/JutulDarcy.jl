function darcy_permeability_potential_differences(
        face,
        state,
        model,
        flux_type,
        mpfa::D,
        upw,
        phases = eachphase(model.system)
    ) where {D<:Jutul.NFVM.NFVMDiscretization}
    pc, ref_index = capillary_pressure(model, state)
    p = state.Pressure
    # if haskey(state, :PermeabilityMultiplier)
    #     K_mul = state.PermeabilityMultiplier
    #     m = face_average(c -> @inbounds K_mul[c], mpfa)
    # end
    m = 1.0
    pot = state.PhasePotentials
    dens = state.PhaseMassDensities
    z = state.AdjustedCellDepths
    # grad(p - rho g z) = grad(p) - grad(rho) g z - rho g grad(z)
    # -> grad(p) - rho g grad(z) = grad(p - rho g z) + grad(rho) g z
    l, r = Jutul.cell_pair(mpfa)
    z_avg = (z[l] + z[r])/2.0
    q = map(ph -> nfvm_potential_difference(pot, dens, z_avg, mpfa, ph, m), phases)
    return q
end

function nfvm_potential_difference(pot, dens, z_avg, mpfa, phase, m = 1.0)
    g = gravity_constant
    ∇pot = Jutul.NFVM.evaluate_flux(pot, mpfa, phase)
    ∇rho = Jutul.NFVM.evaluate_flux(dens, mpfa, phase)
    return -m*(∇pot + ∇rho*g*z_avg)
end

@inline function JutulDarcy.gradient(X::AbstractVector, hf::Jutul.NFVM.NFVMDiscretization)
    l, r = Jutul.cell_pair(hf)
    return @inbounds X[r] - X[l]
end

@inline function JutulDarcy.gradient(X::AbstractMatrix, i, hf::Jutul.NFVM.NFVMDiscretization)
    l, r = Jutul.cell_pair(hf)
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
