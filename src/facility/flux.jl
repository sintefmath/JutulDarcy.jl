function select_secondary_variables_flow_type!(S, domain, system, formulation, flow_type::TotalMassVelocityMassFractionsFlow)
    S[:TotalMass] = TotalMass()
end
include_face_sign(::TotalMassVelocityMassFractionsFlow) = true


"""
Total velocity half face flux with mass fractions for TotalMassVelocityMassFractionsFlow
"""
function update_half_face_flux!(law::ConservationLaw, storage, model, dt, flowd::TwoPointPotentialFlow{U, K, T}) where {U,K,T<:TotalMassVelocityMassFractionsFlow}
    state = storage.state
    masses = state.TotalMasses
    total = state.TotalMass
    v = state.TotalMassFlux

    flux_cells = get_entries(law.half_face_flux_cells)
    flux_faces = get_entries(law.half_face_flux_faces)

    conn_data = law.flow_discretization.conn_data
    N = model.domain.grid.neighborship
    update_fluxes_total_mass_velocity_cells!(flux_cells, conn_data, masses, total, v)
    update_fluxes_total_mass_velocity_faces!(flux_faces, N, masses, total, v)
end

function update_fluxes_total_mass_velocity_cells!(flux, conn_data, masses, total, v)
    function q(c, masses, total, v, phno)
        masses_i = view(masses, phno, :)
        f = c.face
        # sign == 1 if current (self) cell is the first in N[:, i], -1 otherwise.
        # Flux i is positive if going from N[1, i] (L) to N[2, i] (R) and negative 
        # otherwise.
        #
        # This means that we need to flip the sign for the purpose of cell-cell
        # fluxes since their convention is to have negative fluxes going out from
        # the self cell.
        s = -c.face_sign
        vi = s*v[f]
        return half_face_fluxes_total_mass_velocity(c.self, c.other, masses_i, total, vi)
    end
    @tullio flux[phno, i] = q(conn_data[i], masses, total, v, phno)
end

function update_fluxes_total_mass_velocity_faces!(flux, N, masses, total, v)
    @tullio flux[ph, f] = half_face_fluxes_total_mass_velocity_face(N[1, f], N[2, f], view(masses, ph, :), total, v[f])
end

function half_face_fluxes_total_mass_velocity(self, other, masses, total, v)
    if v < 0
        # Flux is leaving the cell
        x = masses[self]/total[self]
    else
        # Flux is entering the cell
        x = value(masses[other])/value(total[other])
    end
    return x*value(v)
end

function half_face_fluxes_total_mass_velocity_face(left, right, masses, total, v)
    # Note the different signs. The above function (for cells) compute the half face flux
    # and recieve the signed flux going into or out of the cell. For the half face velocity
    # we have a single velocity, and the convention is to take the left cell to be upstream 
    # for a positive flux.
    if v > 0
        # Flow from left to right
        x = value(masses[left])/value(total[left])
    else
        # Flow from right to left
        x = value(masses[right])/value(total[right])
    end
    return x*v
end
