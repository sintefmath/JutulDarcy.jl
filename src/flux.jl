export DarcyMassMobilityFlow
struct DarcyMassMobilityFlow <: FlowType end
struct DarcyMassMobilityFlowFused <: FlowType end


function select_secondary_variables_flow_type!(S, domain, system, formulation, flow_type::Union{DarcyMassMobilityFlow, DarcyMassMobilityFlowFused})
    if isa(system, SinglePhaseSystem)
        S[:RelativePermeabilities] = ConstantVariables([1.0])
    else
        S[:RelativePermeabilities] = BrooksCoreyRelPerm(system)
    end
end


function single_unique_potential(model)
    # We should add capillary pressure here ventually
    return model.domain.discretizations.mass_flow.gravity
end

"""
Half face Darcy flux with separate potential.
"""
function update_half_face_flux!(law::ConservationLaw, storage, model, dt, flow_type)
    state = storage.state
    param = storage.parameters
    flux = get_entries(law.half_face_flux_cells)
    update_half_face_flux!(flux, state, model, param, dt, flow_type)
end

function update_half_face_flux!(flux::AbstractArray, state, model, param, dt, flow_disc::TwoPointPotentialFlow{U, K, T}) where {U,K,T<:DarcyMassMobilityFlow}
    rho, kr, mu, p = state.PhaseMassDensities, state.RelativePermeabilities, state.PhaseViscosities, state.Pressure
    pc, ref_index = capillary_pressure(model, state)
    conn_data = flow_disc.conn_data
    ctx = model.context
    single_potential = !flow_disc.gravity && isnothing(pc)
    phases = get_phases(model.system)
    update_immiscible_fluxes!(flux, ctx, conn_data, rho, kr, mu, p, pc, ref_index, phases, single_potential)
end

function update_immiscible_fluxes!(flux, context, conn_data, rho, kr, mu, p, pc, ref_index, phases, single_potential)
    mb = minbatch(context)
    if single_potential
        @batch minbatch = mb for i in eachindex(conn_data)
            immiscible_multiphase_flux_single_pot!(flux, i, conn_data, phases, p, rho, kr, mu)
        end
    else
        # Multiphase potential
        @batch minbatch = mb for i in eachindex(conn_data)
            immiscible_multiphase_flux_multi_pot!(flux, i, conn_data, phases, p, rho, kr, mu, pc, ref_index)
        end
    end
end

function update_immiscible_fluxes!(flux, context::SingleCUDAContext, conn_data, rho, kr, mu, p, pc, ref_index, single_potential)
    if single_potential
        # Scalar potential
        @tullio flux[ph, i] = immiscible_flux_for_phase_single_pot(conn_data[i], ph, p, rho, kr, mu)
    else
        # Multiphase potential
        @tullio flux[ph, i] = immiscible_flux_for_phase_multi_pot(conn_data[i], ph, p, rho, kr, mu, pc, ref_index)
    end
end

function immiscible_multiphase_flux_multi_pot!(flux, conn_i, conn_data, phases, p, rho, kr, mu, pc, ref_index)
    @inbounds cd = conn_data[conn_i]
    self, other, g??z, T = cd.self, cd.other, cd.gdz, cd.T
    ??? = (x) -> local_ad(x, self)

    kr = ???(kr)
    mu = ???(mu)
    rho = ???(rho)
    p = ???(p)
    pc = ???(pc)
    for ph in eachindex(phases)
        @inbounds flux[ph, conn_i] = immiscible_flux_gravity(self, other, ph, kr, mu, rho, p, pc, T, g??z, ref_index)
    end
end

@inline function immiscible_flux_for_phase_multi_pot(cd, ph, p, rho, kr, mu, pc, ref_index)
    c, i, g??z, T = cd.self, cd.other, cd.gdz, cd.T
    ??? = (x) -> local_ad(x, c)
    return immiscible_flux_gravity(c, i, ph, ???(kr), ???(mu), ???(rho), ???(p), ???(pc), T, g??z, ref_index)
end

@inline function immiscible_flux_gravity(c, i, ph, k???, ??, ??, P, pc, T, g??z, ref_index)
    @inbounds ??_c = ??[ph, c]
    @inbounds ??_i = ??[ph, i]
    ??_avg = (??_i + ??_c)/2
    ??pc = capillary_gradient(pc, c, i, ph, ref_index)
    @inbounds ??p = P[c] - P[i] + ??pc
    ?? = -T*(??p + g??z*??_avg)
    if ?? < 0
        # Flux is leaving the cell
        up = c
        ?? = ??_c
    else
        # Flux is entering the cell
        up = i
        ?? = ??_i
    end
    @inbounds ??????? = ??*k???[ph, up]/??[ph, up]
    return ???????*??
end


function immiscible_multiphase_flux_single_pot!(flux, conn_i, conn_data, phases, p, rho, kr, mu)
    @inbounds cd = conn_data[conn_i]
    self, other, T = cd.self, cd.other, cd.T
    ??? = (x) -> local_ad(x, self)

    kr = ???(kr)
    mu = ???(mu)
    rho = ???(rho)
    p = ???(p)
    for ph in eachindex(phases)
        @inbounds flux[ph, conn_i] = immiscible_flux_no_gravity(self, other, ph, kr, mu, rho, p, T)
    end
end

function immiscible_flux_for_phase_single_pot(cd, ph, p, rho, kr, mu)
    c, i, T = cd.self, cd.other, cd.T
    ??? = (x) -> local_ad(x, c)
    return immiscible_flux_no_gravity(c, i, ph, ???(kr), ???(mu), ???(rho), ???(p), T)
end

function immiscible_flux_no_gravity(c, i, ph, k???, ??, ??, P, T)
    @inbounds ?? = -T*(P[c] - P[i])
    if ?? < 0
        # Flux is leaving the cell
        up_c = c
    else
        # Flux is entering the cell
        up_c = i
    end
    @inbounds ??????? = ??[ph, up_c]*k???[ph, up_c]/??[ph, up_c]
    return ???????*??
end


capillary_gradient(::Nothing, c_l, c_r, ph, ph_ref) = 0.0
function capillary_gradient(pc, c_l, c_r, ph, ph_ref)
    if ph == ph_ref
        ??p_c = 0.0
    elseif ph < ph_ref
        ??p_c = pc[ph, c_l] - pc[ph, c_r]
    else
        ??p_c = pc[ph-1, c_l] - pc[ph-1, c_r]
    end
end



"""
TPFA KGrad(p) without gravity. (Outer version, with conn_data input)
"""
@inline function half_face_two_point_kgradp(conn_data::NamedTuple, p::AbstractArray)
    half_face_two_point_kgradp(conn_data.self, conn_data.other, conn_data.T, p)
end

"""
TPFA KGrad(p) without gravity. (Inner version, with explicit inputs)
"""
@inline function half_face_two_point_kgradp(c_self::I, c_other::I, T, p::AbstractArray{R}) where {R<:Real, I<:Integer}
    return -T*(p[c_self] - value(p[c_other]))
end

"""
TPFA-SPU Mobility * KGrad(p) without gravity. (Outer version, with conn_data input)
"""
@inline function half_face_two_point_flux_fused(conn_data, p, ??)
    return half_face_two_point_flux_fused(conn_data.self, conn_data.other, conn_data.T, p, ??)
end

"""
TPFA-SPU Mobility * KGrad(p) without gravity. (Inner version, with explicit inputs)
"""
@inline function half_face_two_point_flux_fused(c_self, c_other, T, p, ??)
    ?? = half_face_two_point_kgradp(c_self, c_other, T, p)
    ????? = spu_upwind(c_self, c_other, ??, ??)
    return ?????*??
end

"""
TPFA-SPU Mobility * (KGrad(p) + G). (Outer version, with conn_data input)
"""
@inline function half_face_two_point_flux_fused_gravity(conn_data, p, ??, density)
    return half_face_two_point_flux_fused_gravity(conn_data.self, conn_data.other, conn_data.T, p, ??, conn_data.gdz, density)
end

"""
TPFA-SPU Mobility * (KGrad(p) + G). (Inner version, with explicit inputs)
"""
@inline function half_face_two_point_flux_fused_gravity(c_self, c_other, T, p, ??, g??z, density)
    ?? = half_face_two_point_kgradp_gravity(c_self, c_other, T, p, g??z, density)
    return spu_upwind_mult(c_self, c_other, ??, ??)
end


"""
Two point Darcy flux with gravity - outer version that takes in NamedTuple for static parameters
"""
@inline function half_face_two_point_kgradp_gravity(conn_data::NamedTuple, p, density)
    return half_face_two_point_kgradp_gravity(conn_data.self, conn_data.other, conn_data.T, p, conn_data.gdz, density)
end

"""
Two point Darcy flux with gravity - inner version that takes in cells and transmissibily explicitly
"""
@inline function half_face_two_point_kgradp_gravity(c_self::I, c_other::I, T, p::AbstractArray{R}, g??z, ??::AbstractArray{R}) where {R<:Real, I<:Integer}
    v = -T*Jutul.two_point_potential_drop_half_face(c_self, c_other, p, g??z, ??)
    return v
end

@inline function phase_mass_flux(??, c, i, ??, kr, ??, ph)
    upc = upwind_cell(??, c, i)
    @inbounds F = ??[ph, upc]*(kr[ph, upc]/??[ph, upc])*??
    return (F, upc)
end

@inline function upwind_cell(pot, l, r)
    if pot < 0
        c = l
    else
        c = r
    end
end

@inline function saturation_averaged_density(??, ph, sat, c1, c2)
    @inbounds ??_1 = ??[ph, c1]
    @inbounds ??_2 = ??[ph, c2]
    @inbounds S_1 = sat[ph, c1]
    @inbounds S_2 = sat[ph, c2]

    avg = (??_1*S_1 + ??_2*S_2)/max(S_1 + S_2, 1e-12)
    return avg
end
