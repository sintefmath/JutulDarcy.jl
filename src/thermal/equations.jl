@inline function Jutul.face_flux!(Q, left, right, face, face_sign, eq::ConservationLaw{:TotalThermalEnergy, <:Any}, state, model, dt, flow_disc::TwoPointPotentialFlowHardCoded)
    # Specific version for tpfa flux
    # TODO: Add general version for thermal
    grad = TPFA(left, right, face_sign)
    upw = SPU(left, right)
    # TODO: This could be inconsistent if we use flux_type for something.
    flux_type = Jutul.flux_type(eq)
    q = thermal_heat_flux(face, state, model, grad, upw, flux_type)
    return setindex(Q, q, 1)
end

@inline function Jutul.face_flux!(q_i, face, eq::ConservationLaw{:TotalThermalEnergy, <:Any}, state, model, dt, flow_disc::PotentialFlow, ldisc)
    # Inner version, for generic flux
    kgrad, upw = ldisc.face_disc(face)
    ft = Jutul.flux_type(eq)
    q = thermal_heat_flux(face, state, model, kgrad, upw, ft)
    return setindex(q_i, q, 1)
end

"""
    thermal_heat_flux(face, state, model, grad, upw, flux_type)

Calculate the thermal heat flux for a given face in a thermal model.

# Arguments
- `face`: The face for which the heat flux is being calculated.
- `state`: The current state of the system.
- `model`: The thermal model being used.
- `grad`: The gradient operator.
- `upw`: Upwind scheme operator.
- `flux_type`: The type of flux calculation to be used.

# Returns
- The calculated thermal heat flux for the given face.
"""
function thermal_heat_flux(face, state, model, grad, upw, flux_type)
    T = state.Temperature
    H_f = state.FluidEnthalpy
    θ_r = state.RockThermalTransmissibilites
    S = state.Saturations
    nph = number_of_phases(model.system)

    convective_flux = 0.0
    mass_fluxes = darcy_phase_mass_fluxes(face, state, model, flux_type, grad, upw)
    for α in 1:nph
        F_α = mass_fluxes[α]
        H_face_α = phase_upwind(upw, H_f, α, F_α)
        convective_flux += H_face_α*F_α
    end

    # Compute fluid thermal transmissibilities dynamically when conductivity is available.
    has_dynamic_λf = haskey(model.secondary_variables, :FluidThermalConductivity) && haskey(state, :FluidThermalConductivity)
    if has_dynamic_λf
        λ_f = state.FluidThermalConductivity
        ϕ = state.FluidVolume./state.BulkVolume
        domain = model.data_domain
        dim = size(domain[:cell_centroids], 1)
    else
        θ_f = state.FluidThermalTransmissibilites
    end

    θ = θ_r[face]
    for α in 1:nph
        if has_dynamic_λf
            grad isa TPFA || throw(ArgumentError("TPFA gradient expected for dynamic fluid thermal conductivity"))
            
            left, right = domain[:neighbors][:, face]
            # left == grad.left || throw(ArgumentError("Left cell mismatch in TPFA gradient"))
            # right == grad.right || throw(ArgumentError("Right cell mismatch in TPFA gradient"))

            A = domain[:areas][face]
            N = domain[:normals][:, face]

            den = 0.0
            for side in (left, right)
                C = domain[:face_centroids][:, face] - domain[:cell_centroids][:, side]
                sgn = ifelse(side == left, 1.0, -1.0)
                λ = ϕ[side].*fluid_phase_cell_value(λ_f, α, side)
                λ = Jutul.expand_perm(λ, Val(dim))
                θ_hf = Jutul.half_face_trans(A, λ, C, sgn*N)
                den += 1/θ_hf
            end
            θ_f = 1/den

            # λ_face_α = @inbounds λ_f[α, face]
        else
            θ_f = θ_f[α, face]
        end
        θ += θ_f*phase_face_average(S, grad, α)
    end

    if haskey(state, :ThermalDispersivity)
        neighbors = model.data_domain[:neighbors]
        # q_phases = face_flux_helper(face, neighbors, state, model, false)
        # qT = sum(q_phases)
        domain = model.data_domain
        dim = size(domain[:cell_centroids], 1)
        left, right = neighbors[:, face]
        A = domain[:areas][face]
        Nf = domain[:normals][:, face]
        den = 0.0
        Id = SMatrix{3,3,Float64}(I)
        for side in (left, right)
            # Compute the dispersive flux contribution for each cell adjacent to the face
            v = velocity_from_state(side, state, model; is_mass = false)
            v = SVector{dim}(v)
            v² = sqrt(v'*v)
            if v² == 0.0
                continue
            end
            α_L = state.ThermalDispersivity[1, side]
            α_T = state.ThermalDispersivity[2, side]
            ρ_f = state.PhaseMassDensities[1, side]
            Cp_f = domain[:component_heat_capacity][side]

            # D = @SVector [
            #     α_L*v[1].^2/v_norm + α_T*(v[2].^2 + v[3].^2)/v_norm,
            #     α_L*v[2].^2/v_norm + α_T*(v[1].^2 + v[3].^2)/v_norm,
            #     α_L*v[3].^2/v_norm + α_T*(v[1].^2 + v[2].^2)/v_norm
            # ]
            ρ_f = state.PhaseMassDensities[1, side]
            Cp_f = domain[:component_heat_capacity][side]

            D = ρ_f*Cp_f*(α_T*Id + (α_L - α_T)*(v*v')/(v'*v)).*sqrt((v'*v))
            # D .*= ρ_f*Cp_f
            # D = ρ*Cp*(αT*Id + (αL - αT)*(v*v')/(v'*v)).*sqrt((v'*v))
            C = domain[:face_centroids][:, face] - domain[:cell_centroids][:, side]
            sgn = ifelse(side == left, 1.0, -1.0)
            # D = Jutul.expand_perm(D, Val(dim))
            # display(value(D))
            θ_hf = Jutul.half_face_trans(A, D, C, sgn*Nf)
            den += 1/θ_hf
        end

        if den == 0.0
            θ_d = 0.0
        else
            θ_d = 1/den
        end
        # println("Convective contribution: ", value(den))
        # println("Dispersive contribution: ", value(θ_d))
        θ += θ_d

    end

    conductive_flux = -θ*gradient(T, grad)
    return conductive_flux + convective_flux
end

function fluid_phase_cell_value(λ_f, α, cell)
    if λ_f isa AbstractVector
        return @inbounds λ_f[cell]
    elseif λ_f isa AbstractArray
        return @inbounds λ_f[α, cell]
    end
end

function face_flux_helper(face, N, state, model, is_mass::Bool)
    l = N[1, face]
    r = N[2, face]
    # TODO: This assumes the default discretizations.
    # This should be generalized to allow for different flux types.
    tpfa = TPFA(l, r, 1)
    upw = SPU(l, r)
    f_t = Jutul.DefaultFlux()
    v_face = JutulDarcy.darcy_phase_volume_fluxes(face, state, model, f_t, tpfa, upw)
    return v_face
end

function velocity_from_state(cell, state, model; is_mass::Bool = false)
    if haskey(state, :TotalDarcyVelocity)
        return @view state.TotalDarcyVelocity[:, cell]
    else
        return 0.0#velocity_from_face_fluxes(cell, state, model; is_mass = is_mass)
    end
end

function velocity_from_face_fluxes(cell, state, model; is_mass::Bool = false)
    domain = model.data_domain
    if cell <= 0
        dim = size(domain[:cell_centroids], 1)
        return zeros(dim)
    end

    neighbors = domain[:neighbors]
    nc = size(domain[:cell_centroids], 2)
    faces, facepos = get_facepos(neighbors, nc)
    facesigns = Jutul.get_facesigns(neighbors, faces, facepos, nc)

    cc = domain[:cell_centroids][:, cell]
    v = nothing
    for fpos = facepos[cell]:(facepos[cell+1]-1)
        face = faces[fpos]
        sgn = facesigns[fpos]
        fc = domain[:face_centroids][:, face]
        q = sum(face_flux_helper(face, neighbors, state, model, is_mass))*sgn
        dq = q.*(fc .- cc)
        if isnothing(v)
            v = copy(dq)
        else
            v .+= dq
        end
    end
    if isnothing(v)
        return zeros(eltype(cc), length(cc))
    end
    v ./= domain[:volumes][cell]
    return v
end

"""
    Jutul.convergence_criterion(model, storage, eq::ConservationLaw{:TotalThermalEnergy}, eq_s, r; dt = 1.0, update_report = missing)

Calculate the convergence criterion for the total thermal energy conservation law.

# Arguments
- `model`: The model object containing the simulation parameters and state.
- `storage`: The storage object used to keep track of intermediate results.
- `eq::ConservationLaw{:TotalThermalEnergy}`: The conservation law for total thermal energy.
- `eq_s`: The storage of the conservation law equation.
- `r`: The residual of the conservation law equation.
- `dt`: The time step size in seconds (default is 1.0).
- `update_report`: An optional argument for updating the report (default is `missing`).

# Returns
- The convergence criterion values for the total thermal energy conservation law (maximum).
"""
function Jutul.convergence_criterion(model, storage, eq::ConservationLaw{:TotalThermalEnergy}, eq_s, r; dt = 1.0, update_report = missing)
    a = active_entities(model.domain, Cells())
    E0 = storage.state0.TotalThermalEnergy
    eb, cnv, Etot = 0.0, -Inf, 0.0
    ΔT = temperature_increment(model, storage.state, update_report)
    for (i, c) in enumerate(a)
        eb += r[i]
        cnv = max(cnv, abs(r[i])*dt/value(E0[c]))
        Etot += value(E0[c])
    end
    eb = abs(eb)*dt/Etot

    return (
        CNV = (errors = (cnv, ), names = ("Max", )),
        EB = (errors = (eb, ), names = ("Energy balance", )), 
        increment_dT = (errors = (ΔT, ), names = ("ΔT", )),
        )
end

function temperature_increment(model, state, update_report)
    t_report = update_report[:Temperature]
    if haskey(t_report, :max)
        v = t_report.max
    else
        v = 1.0
    end
    return v
end

function temperature_increment(model, state, update_report::Missing)
    return 1.0
end
