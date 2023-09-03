function flash_wellstream_at_surface(well_model, sys::S, well_state, rhoS, cond = default_surface_cond()) where S<:MultiComponentSystem
    total_masses = well_state.TotalMasses
    T = eltype(total_masses)
    nph = length(rhoS)
    rho = zeros(T, nph)
    volfrac = zeros(T, nph)
    eos = sys.equation_of_state
    nc = MultiComponentFlash.number_of_components(eos)

    if nph == 3
        a, l, v = phase_indices(sys)
        rhoWS = rhoS[a]
        # Convert to surface conditions
        rhoWW = well_state.PhaseMassDensities[a, 1]
        S_other = well_state.Saturations[a, 1]*rhoWW/rhoWS
        # Surface density given for aqueous phase
        rho[a] = rhoWS
        volfrac[a] = S_other
        n = nc + 1
    else
        l, v = phase_indices(sys)
        n = nc
        S_other = zero(T)
    end
    buf = InPlaceFlashBuffer(nc)

    Pressure = cond.p
    Temperature = cond.T

    z = SVector{nc}(well_state.OverallMoleFractions[:, 1])
    m = SSIFlash()
    fstorage = flash_storage(eos, method = m, inc_jac = true, diff_externals = true, npartials = n, static_size = true)
    update_flash_buffer!(buf, eos, Pressure, Temperature, z)

    f = FlashedMixture2Phase(eos, T)
    x = f.liquid.mole_fractions
    y = f.vapor.mole_fractions
    forces = buf.forces

    result = update_flash_result(fstorage, m, eos, f.K, x, y, buf.z, forces, Pressure, Temperature, z)

    rho[l], rho[v] = mass_densities(eos, Pressure, Temperature, result)
    rem = one(T) - S_other
    S_l, S_v = phase_saturations(eos, Pressure, Temperature, result)

    volfrac[l] = rem*S_l
    volfrac[v] = rem*S_v
    rho = tuple(rho...)
    return (rho, volfrac)
end

Base.@propagate_inbounds function multisegment_well_perforation_flux!(out, sys::CompositionalSystem, state_res, state_well, rhoS, conn)
    rc = conn.reservoir
    wc = conn.well

    μ = state_res.PhaseViscosities
    kr = state_res.RelativePermeabilities
    ρ = state_res.PhaseMassDensities
    X = state_res.LiquidMassFractions
    Y = state_res.VaporMassFractions

    ρ_w = state_well.PhaseMassDensities
    s_w = state_well.Saturations
    X_w = state_well.LiquidMassFractions
    Y_w = state_well.VaporMassFractions

    nc = size(X, 1)
    nph = size(μ, 1)
    has_water = nph == 3
    phase_ix = phase_indices(sys)

    if has_water
        A, L, V = phase_ix
    else
        L, V = phase_ix
    end

    mob(ph) = kr[ph, rc]/μ[ph, rc]
    λ_t = zero(eltype(kr))
    for ph in 1:nph
        λ_t += mob(ph)
    end

    function phase_mass_flux(ph)
        dp = perforation_phase_potential_difference(conn, state_res, state_well, ph)
        # dp is pressure difference from reservoir to well. If it is negative,
        # we are injecting into the reservoir.
        if dp < 0
            mob_ph = λ_t*s_w[ph, wc]
            dens_ph = ρ_w[ph, wc]
        else
            mob_ph = mob(ph)
            dens_ph = ρ[ph, rc]
        end
        Q = mob_ph*dens_ph*dp
        return (Q, dp)
    end

    Q_l, dp_l = phase_mass_flux(L)
    Q_v, dp_v = phase_mass_flux(V)

    @inbounds for c in 1:nc
        if dp_l < 0
            # Injection
            X_upw = X_w[c, wc]
        else
            X_upw = X[c, rc]
        end
        if dp_v < 0
            Y_upw = Y_w[c, wc]
        else
            Y_upw = Y[c, rc]
        end
        out[c] = Q_l*X_upw + Q_v*Y_upw
    end
    if has_water
        out[nc+1], _ = phase_mass_flux(A)
    end
end


function separator_surface_flash!(var, model, system::MultiPhaseCompositionalSystemLV, state)
    # For each stage we need to keep track of:
    # Total mass / moles
    # Composition
    # Can then flash w.r.t. that set of conditions
    eos = system.equation_of_state
    nc = MultiComponentFlash.number_of_components(eos)
    z = SVector{nc}(state.OverallMoleFractions[:, 1])
    F = get_separator_flash_buffer(var, system, z)

    streams, moles = get_separator_buffers(var, z)
    moles[1] = z
    streams[1] = z
    n_stage = length(var.separator_conditions)
    for i in 1:n_stage
        cond = var.separator_conditions[i]

        liquid_stream, vapor_stream = flash_stream!(streams[i], moles[i], F, eos, cond)
        if i == n_stage
            # Quit here.
        end
    end
    @info "??" z
    error("Separator system found, not implemented yet")
end

function flash_stream!(stream, moles::SVector{N, T}, flash, eos, cond) where {N, T}
    fstorage, buf, f = flash

    z = moles./sum(moles)
    Pressure = cond.p
    Temperature = cond.T
    update_flash_buffer!(buf, eos, Pressure, Temperature, z)

    x = f.liquid.mole_fractions
    y = f.vapor.mole_fractions
    forces = buf.forces

    result = update_flash_result(fstorage, SSIFlash(), eos, f.K, x, y, buf.z, forces, Pressure, Temperature, z)

    rho, rho = mass_densities(eos, Pressure, Temperature, result)
    S_l, S_v = phase_saturations(eos, Pressure, Temperature, result)

    return (q_l, q_v)
end

function get_separator_flash_buffer(var, system, z)
    s = var.storage
    if !haskey(s, :flash_storage) || true
        nc = number_of_components(system)
        n = nc + has_other_phase(system)
        eos = system.equation_of_state
        m = SSIFlash()
        buf = InPlaceFlashBuffer(nc)
        T = eltype(z)
        f = FlashedMixture2Phase(eos, T)
        fstorage = flash_storage(eos, method = m, inc_jac = true, diff_externals = true, npartials = n, static_size = true)    
        s[:flash_storage] = (fstorage, buf, f)
    end
    return s[:flash_storage]
end

function get_separator_buffers(var, z::T) where T
    # TODO: Deal with varying AD inputs in this function and the flash one.
    s = var.storage
    n = length(var.separator_conditions)
    if !haskey(s, :mole_stages) || true
        s[:mole_stages] = zeros(T, n)
    end
    if !haskey(s, :mole_streams) || true
        s[:mole_streams] = zeros(T, n)
    end
    moles = s[:mole_stages]
    @. moles = zero(T)
    streams = s[:mole_streams]
    @. streams = zero(T)
    return (streams, moles)
end
