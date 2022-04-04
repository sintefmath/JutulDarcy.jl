function apply_forces_to_equation!(storage, model::SimulationModel{D, S}, eq::ConservationLaw, force::V, time) where {V <: AbstractVector{SourceTerm{I, F, T}}, D, S<:CompositionalSystem} where {I, F, T}
    acc = get_diagonal_entries(eq)
    state = storage.state
    kr = state.RelativePermeabilities
    rho = state.PhaseMassDensities
    mu = state.PhaseViscosities
    X = state.LiquidMassFractions
    Y = state.VaporMassFractions
    Sat = state.Saturations
    FR = state.FlashResults
    rhoS = get_reference_densities(model, storage)
    insert_component_sources!(acc, Sat, kr, mu, FR, X, Y, rho, rhoS, force)
end

function insert_component_sources!(acc, S, kr, mu, F, X, Y, rho, rhoS, sources)
    ncomp = size(acc, 1)
    for src in sources
        for c = 1:ncomp
            @inbounds acc[c, src.cell] -= component_source(src, S, kr, mu, F, X, Y, rho, rhoS, c)
        end
    end
end

function component_source(src, S, kr, mu, F, X, Y, rho, rhoS, c)
    # Treat inflow as volumetric with respect to surface conditions
    # Treat outflow as volumetric sources
    v = src.value
    cell = src.cell
    # MassSource -> source terms are already masses
    # VolumeSource -> source terms are volumes at medium conditions
    t = src.type
    f_c = src.fractional_flow[c]
    if v > 0
        if t == MassSource
            # Source term is already as masses. Injection is straightforward, 
            # production is weighted according to mobility*density
            f = f_c
        elseif t == StandardVolumeSource
            f = rhoS[1]*f_c
        elseif t == VolumeSource
            # Saturation-averaged density
            ρ_mix = S[1, cell]*rho[1, cell] + S[2, cell]*rho[2, cell]
            f = f_c*ρ_mix
        else
            error("Not supported source term type: $t")
        end
    else
        f = compositional_out_f(kr, mu, X, Y, rho, rhoS, c, cell, t, F[cell].state)
    end
    q = v*f
    return q
end

function compositional_out_f(kr, mu, X, Y, rho, rhoS, c, cell, t, state)
    λ_l = local_mobility(kr, mu, 1, cell)
    λ_v = local_mobility(kr, mu, 2, cell)

    ρ_l = rho[1, cell]
    ρ_v = rho[2, cell]

    x = X[c, cell]
    y = Y[c, cell]

    if t == MassSource
        m_l = ρ_l*λ_l
        m_v = ρ_v*λ_v
        m_t = m_l + m_v
        f = (m_l*x + m_v*y)/m_t
    elseif t == VolumeSource
        λ_t = λ_l + λ_v
        f = (λ_l*x*ρ_l + λ_v*y*ρ_v)/λ_t
    elseif t == StandardVolumeSource
        ρ_ls = rhoS[1]
        ρ_vs = rhoS[2]
    
        λ_t = λ_l + λ_v
        f_l = λ_l/λ_t
        f_v = λ_v/λ_t
        f = ρ_ls*x*f_l + ρ_vs*y*f_v
    end
    return f
end
