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

function Jutul.select_equations!(eqs, ::TransportFormulation, model::TransportModel)
    lbl = :mass_conservation
    @assert haskey(eqs, lbl)
    eq = eqs[lbl]

    s = conserved_symbol(eq)
    @assert s == :TotalMasses
    N = number_of_equations_per_entity(model, eq)
    eqs[lbl] = ConservationLaw(eq.flow_discretization, s, N, flux = TotalSaturationFlux())
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

    left_mob = map(phase -> state.PhaseMobilities[phase, l], phases)
    right_mob = map(phase -> state.PhaseMobilities[phase, r], phases)

    T_f = effective_transmissibility(state, face, kgrad)
    gΔz = effective_gravity_difference(state, face, kgrad)
    pc, ref_index = capillary_pressure(model, state)

    @inline function bouyancy_and_capillary_term(phase)
        Δpc = capillary_gradient(pc, kgrad, phase, ref_index)
        ρ_avg = face_average_density(model, state, kgrad, phase)
        return -(gΔz*ρ_avg + Δpc)
        # return -T_f*(∇p + Δpc + gΔz*ρ_avg)
    end

    G = map(bouyancy_and_capillary_term, phases)

    if false
        @assert number_of_phases(model.system) == 2 "This debug option is only implemented for two phases"
        mob_1 = JutulDarcy.phase_upwind(upw, state.PhaseMobilities, 1, V_t)
        mob_2 = JutulDarcy.phase_upwind(upw, state.PhaseMobilities, 2, V_t)
        mob_t = mob_1 + mob_2
        return (1/mob_t*V_t, 1/mob_t*V_t)
    end
    return phase_potential_upwind_potential_differences(V_t, T_f, G, left_mob, right_mob)
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
        q_p *= state_s.TotalSaturation[conn.well]
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

