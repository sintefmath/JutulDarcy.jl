function Jutul.select_primary_variables!(vars, ::TransportFormulation, model::TransportModel)
    delete!(vars, :Pressure)
    vars[:TotalSaturation] = TotalSaturation()
end

function Jutul.select_parameters!(vars, ::TransportFormulation, model::TransportModel)
    vars[:Pressure] = Pressure()
    if JutulDarcy.model_or_domain_is_well(model)
        vars[:PerforationTotalVolumetricFlux] = PerforationTotalVolumetricFlux()
    else
        vars[:TotalVolumetricFlux] = TotalVolumetricFlux()
    end
end

function Jutul.select_equations!(eqs, f::TransportFormulation, model::TransportModel)
    lbl = :mass_conservation
    @assert haskey(eqs, lbl)
    eq = eqs[lbl]

    s = conserved_symbol(eq)
    @assert s == :TotalMasses
    N = number_of_equations_per_entity(model, eq)
    eqs[lbl] = ConservationLaw(eq.flow_discretization, s, N, flux = TotalSaturationFlux(f))
end

@inline function JutulDarcy.darcy_permeability_potential_differences(
        face,
        state,
        model,
        flux_type::TotalSaturationFlux,
        kgrad::TPFA,
        upw,
        phases = JutulDarcy.eachphase(model.system)
    )
    V_t = kgrad.face_sign*state.TotalVolumetricFlux[face]
    N = length(phases)
    @assert N == number_of_phases(model.system)
    l, r = Jutul.cell_pair(upw)

    mob(phase, c) = @inbounds state.PhaseMobilities[phase, c]
    left_mob = map(phase -> mob(phase, l), phases)
    right_mob = map(phase -> mob(phase, r), phases)

    T_f = effective_transmissibility(state, face, kgrad)
    gΔz = effective_gravity_difference(state, face, kgrad)
    pc, ref_index = capillary_pressure(model, state)

    bouyancy_term(phase) = -gΔz*face_average_density(model, state, kgrad, phase)
    capillary_term(phase) = -capillary_gradient(pc, kgrad, phase, ref_index)
    @inline function bouyancy_and_capillary_term(phase)
        # Note: These have their signs flipped already!
        Δpc = capillary_term(phase)
        gΔz_ρ_avg = bouyancy_term(phase)
        # Reference:  -T_f*(∇p + Δpc + gΔz*ρ_avg) for fully coupled flux
        return (gΔz_ρ_avg + Δpc)
    end

    scheme = transport_scheme(flux_type)
    if scheme == :ppu
        G = map(bouyancy_and_capillary_term, phases)
        out = phase_potential_upwind_potential_differences(V_t, T_f, G, left_mob, right_mob)
    else
        G_den = map(bouyancy_term, phases)
        G_cap = map(capillary_term, phases)
        pot_den = phase_potential_upwind_potential_differences(V_t, T_f, G_den, left_mob, right_mob)

        if scheme == :ppu_nopc || all(G_cap .== 0.0)
            out = pot_den
        else
            pot_pc = phase_potential_upwind_potential_differences(zero(V_t), T_f, G_cap, left_mob, right_mob)
            out = MultiPotential(pot_den, pot_pc)
        end
    end
    return out
end

function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t,
        state_s, state0_s, 
        model_t::TransportModel, model_s::TransportModel,
        ct::ReservoirFromWellFlowCT, eq, dt, ldisc = local_discretization(ct, i)
    )
    sys = model_t.system
    rhoS = reference_densities(sys)
    conn = cross_term_perforation_get_conn(ct, i, state_s, state_t)
    q_p = state_s[:PerforationTotalVolumetricFlux][i]
    # Call smaller interface that is easy to specialize
    @assert !haskey(state_s, :MassFractions)
    # @assert abs(conn.gdz) < 1e-10 "connection gravity difference should be zero for transport model, was $(conn.gdz)"
    rc = conn.reservoir
    mob = state_t.PhaseMobilities
    mobt = zero(eltype(mob))
    for ph in axes(mob, 1)
        mobt += mob[ph, rc]
    end
    injecting = q_p > 0.0
    if injecting
        # q_p *= state_s.TotalSaturation[conn.well]
    end
    conn = (
        dp = q_p/(mobt*conn.WI),
        WI = conn.WI,
        gdz = 0*conn.gdz,
        well = conn.well,
        perforation = conn.perforation,
        reservoir = conn.reservoir
    )
    JutulDarcy.multisegment_well_perforation_flux!(out, sys, state_t, state_s, rhoS, conn)
    # Add in a hack to ensure that sparsity gets properly detected.
    dp = state_t.TotalSaturation[conn.reservoir] - state_s.TotalSaturation[conn.well]
    out[1] += 0.0*dp
    return out
end

function JutulDarcy.reservoir_linsolve(model::TransportModel, pname = :ilu0;
        rtol = 1e-3,
        solver = :bicgstab,
        kwarg...
    )
    if pname == :ilu0
        prec = Jutul.ILUZeroPreconditioner()
        lsolve = GenericKrylov(solver; preconditioner = prec, rtol = rtol, kwarg...)
    else
        error("Preconditioner $pname not supported for transport model")
    end
    return lsolve
end

