function flash_wellstream_at_surface(well_model::SimulationModel{D, S}, well_state, rhoS) where {D, S<:MultiComponentSystem}
    total_masses = well_state.TotalMasses
    T = eltype(total_masses)
    nph = length(rhoS)
    rho = zeros(T, nph)
    volfrac = zeros(T, nph)
    sys = well_model.system
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

    sc = well_model.domain.grid.surface
    Pressure = sc.p
    Temperature = sc.T

    z = SVector{nc}(well_state.OverallMoleFractions[:, 1])
    m = SSIFlash()
    fstorage = flash_storage(eos, method = m, inc_jac = true, diff_externals = true, npartials = n, static_size = true)
    update_flash_buffer!(buf, eos, Pressure, Temperature, z)

    f = FlashedMixture2Phase(eos, T)
    x = f.liquid.mole_fractions
    y = f.vapor.mole_fractions
    forces = buf.forces

    result = update_flash_result(fstorage, m, buf, eos, f.K, x, y, buf.z, forces, Pressure, Temperature, z)

    rho[l], rho[v] = mass_densities(eos, Pressure, Temperature, result)
    rem = one(T) - S_other
    S_l, S_v = phase_saturations(eos, Pressure, Temperature, result)

    volfrac[l] = rem*S_l
    volfrac[v] = rem*S_v
    rho = tuple(rho...)
    return (rho, volfrac)
end

function flash_wellstream_at_surface(well_model::SimulationModel{D, S}, well_state, rhoS) where {D, S<:ImmiscibleSystem, T}
    tot_mass = well_state.TotalMasses
    vol = map((i, rho) -> tot_mass[i, 1], 1:length(rhoS), rhoS)
    volfrac = vol./sum(vol)
    return (rhoS, volfrac)
end
