@enum PresentPhasesBlackOil OilOnly GasOnly OilAndGas


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
    S[:PhaseMassDensities] = DeckPhaseMassDensities(pvt)
    S[:ShrinkageFactors] = DeckShrinkageFactors(pvt)
    g = physical_representation(model.domain)
    if !(g isa WellDomain)
        S[:SurfaceVolumeMobilities] = SurfaceVolumeMobilities()
    end
    S[:PhaseViscosities] = DeckPhaseViscosities(pvt)
    pvt_regions = reservoir_regions(model, :pvtnum)
    if has_disgas(system)
        S[:Rs] = Rs(system.rs_max, regions = pvt_regions)
    end
    if has_vapoil(system)
        S[:Rv] = Rv(system.rv_max, regions = pvt_regions)
    end
end

get_phases(sys::StandardBlackOilSystem) = sys.phases
number_of_components(sys::StandardBlackOilSystem) = length(get_phases(sys))
phase_indices(sys::StandardBlackOilSystem) = sys.phase_indices

function component_names(sys::StandardBlackOilSystem)
    return phase_names(sys)
end

has_vapoil(::Any) = false
has_disgas(::Any) = false

has_vapoil(::StandardBlackOilSystem) = true
has_disgas(::StandardBlackOilSystem) = true

has_vapoil(::DisgasBlackOilSystem) = false
has_disgas(::VapoilBlackOilSystem) = false

function convergence_criterion(model::SimulationModel{D, S}, storage, eq::ConservationLaw{:TotalMasses}, eq_s, r; dt = 1.0, update_report = missing) where {D, S<:StandardBlackOilSystem}
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

function handle_alternate_primary_variable_spec!(init, found, rmodel, sys::StandardBlackOilSystem)
    # Internal utility to handle non-trivial specification of primary variables
    nph = number_of_phases(sys)
    haskey(init, :Pressure) || error("Primary variable :Pressure is missing from the initial state.")
    pressure = init[:Pressure]
    nc = length(pressure)
    # Check required inputs to define the "black oil unknown"
    has_sat = haskey(init, :Saturations)
    has_box = haskey(init, :BlackOilUnknown)
    has_sosg = haskey(init, :LiquidSaturation) && haskey(init, :VaporSaturation)
    has_sat || has_box || has_sosg || error("Primary variable :Saturations, :LiquidSaturation+:VaporSaturation or :BlackOilUnknown is missing from the initial state.")

    if nph == 3 && !haskey(init, :ImmiscibleSaturation)
        S = init[:Saturations]
        a, l, v = phase_indices(sys)
        sw = S[a, :]
        init[:ImmiscibleSaturation] = sw
        push!(found, :ImmiscibleSaturation)
    end

    if nph == 2
        sw = zeros(nc)
        l, v = phase_indices(sys)
    else
        a, l, v = phase_indices(sys)
        if haskey(init, :ImmiscibleSaturation)
            sw = init[:ImmiscibleSaturation]
        elseif haskey(init, :Saturations)
            sw = init[:Saturations][a, :]
        end
    end

    if !haskey(init, :BlackOilUnknown)
        if has_sosg
            so = init[:LiquidSaturation]
            sg = init[:VaporSaturation]
        else
            so = init[:Saturations][l, :]
            sg = init[:Saturations][v, :]
        end
        F_rs = sys.rs_max
        F_rv = sys.rv_max
        T = promote_type(eltype(so), eltype(sg), eltype(sw), eltype(pressure))
        if has_disgas(sys)
            rs = init[:Rs]
            T = promote_type(eltype(rs), T)
            rs_var = rmodel[:Rs]
        else
            rs = zeros(nc)
            rs_var = nothing
        end
        if has_vapoil(sys)
            rv = init[:Rv]
            T = promote_type(eltype(rv), T)
            rv_var = rmodel[:Rv]
        else
            rv = zeros(nc)
            rv_var = nothing
        end
        so = @. 1.0 - sw - sg
        bo = Vector{BlackOilX{T}}()
        sizehint!(bo, nc)
        for i in 1:nc
            reg_rs = region(rs_var, i)
            reg_rv = region(rv_var, i)
            F_rs_i = table_by_region(F_rs, reg_rs)
            F_rv_i = table_by_region(F_rv, reg_rv)
            sw_i, so_i, sg_i, rs_i, rv_i, pressure_i = promote(sw[i], so[i], sg[i], rs[i], rv[i], pressure[i])
            v = blackoil_unknown_init(F_rs_i, F_rv_i, sw_i, so_i, sg_i, rs_i, rv_i, pressure_i)
            push!(bo, v)
        end
        init[:BlackOilUnknown] = bo
        push!(found, :BlackOilUnknown)
    end
    return init
end
