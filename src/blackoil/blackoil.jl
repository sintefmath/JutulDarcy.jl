@enum PresentPhasesBlackOil OilOnly GasOnly OilAndGas MostlyWater


include("variables/variables.jl")
include("flux.jl")
include("wells.jl")
include("data.jl")
include("utils.jl")

blackoil_formulation(::StandardBlackOilSystem{V, D, W, R, F}) where {V, D, W, R, F} = F

function select_primary_variables!(S, system::BlackOilSystem, model)
    S[:Pressure] = Pressure()
    if has_other_phase(system)
        S[:ImmiscibleSaturation] = ImmiscibleSaturation(ds_max = 0.2)
    end
    bf = blackoil_formulation(system)
    if bf == :varswitch
        S[:BlackOilUnknown] = BlackOilUnknown()
    elseif bf == :zg
        S[:GasMassFraction] = GasMassFraction(dz_max = 0.1)
    else
        error("Unsupported formulation $bf.")
    end
end

function select_secondary_variables!(S, system::BlackOilSystem, model)
    select_default_darcy_secondary_variables!(S, model.domain, model.system, model.formulation)
    S[:Saturations] = Saturations()
    S[:PhaseState] = BlackOilPhaseState()
    spe1_data = blackoil_bench_pvt(:spe1)
    pvt = spe1_data[:pvt]
    S[:PhaseMassDensities] = DeckDensity(pvt)
    S[:ShrinkageFactors] = DeckShrinkageFactors(pvt)
    if !(model.domain.grid isa WellGrid)
        S[:SurfaceVolumeMobilities] = SurfaceVolumeMobilities()
    end
    S[:PhaseViscosities] = DeckViscosity(pvt)
    if has_disgas(system)
        S[:Rs] = Rs()
    end
    if has_vapoil(system)
        S[:Rv] = Rv()
    end
end

get_phases(sys::StandardBlackOilSystem) = sys.phases
number_of_components(sys::StandardBlackOilSystem) = length(get_phases(sys))
phase_indices(sys::StandardBlackOilSystem) = sys.phase_indices

has_vapoil(::Any) = false
has_disgas(::Any) = false

has_vapoil(::StandardBlackOilSystem) = true
has_disgas(::StandardBlackOilSystem) = true

has_vapoil(::DisgasBlackOilSystem) = false
has_disgas(::VapoilBlackOilSystem) = false

function convergence_criterion(model::SimulationModel{D, S}, storage, eq::ConservationLaw, eq_s, r; dt = 1) where {D, S<:StandardBlackOilSystem}
    M = global_map(model.domain)
    v = x -> as_value(Jutul.active_view(x, M, for_variables = false))
    Φ = v(storage.state.FluidVolume)
    b = v(storage.state.ShrinkageFactors)

    sys = model.system
    nph = number_of_phases(sys)
    rhoS = reference_densities(sys)
    cnv, mb = cnv_mb_errors_bo(r, Φ, b, dt, rhoS, Val(nph))

    names = phase_names(model.system)
    R = (CNV = (errors = cnv, names = names),
         MB = (errors = mb, names = names))
    return R
end

function cnv_mb_errors_bo(r, Φ, b, dt, rhoS, ::Val{N}) where N
    nc = length(Φ)
    mb = @MVector zeros(N)
    cnv = @MVector zeros(N)
    avg_B = @MVector zeros(N)

    pv_t = 0.0
    @inbounds for c in 1:nc
        pv_c = Φ[c]
        pv_t += pv_c
        @inbounds for ph = 1:N
            r_ph = r[ph, c]
            b_ph = b[ph, c]
            # MB
            mb[ph] += r_ph
            avg_B[ph] += 1/b_ph
            # CNV
            cnv[ph] = max(cnv[ph], abs(r_ph)/pv_c)
        end
    end
    @inbounds for ph = 1:N
        B = avg_B[ph]/nc
        scale = B*dt/rhoS[ph]
        mb[ph] = scale*abs(mb[ph])/pv_t
        cnv[ph] = scale*abs(cnv[ph])
    end
    return (Tuple(cnv), Tuple(mb))
end
