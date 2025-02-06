import JutulDarcy: region, table_by_region, model_or_domain_is_well, phase_indices, reference_densities, upwind

struct PolymerTracer <: AbstractTracer
    ix::Int
end

function PolymerTracer(system::MultiPhaseSystem, dest_phase = AqueousPhase)
    for (i, ph) in enumerate(get_phases(system))
        if ph isa dest_phase
            return PolymerTracer(i)
        end
    end
    error("Phase $dest_phase not found in system")
end

tracer_phase_indices(t::PolymerTracer) = (t.ix, )

function tracer_total_mass_outer(tracer::PolymerTracer, model, state, concentration, resident_mass_density, pore_vol, cell, index)
    if resident_mass_density <= TRACER_TOL
        v = pore_vol*Jutul.replace_value(concentration, 0.0)
    else
        v = tracer_total_mass(tracer, model, state, concentration, resident_mass_density, pore_vol, cell, index)
    end
    if model_or_domain_is_well(model)
        mass = v
    else
        sys = model.system
        phase_ix = phase_indices(sys)
        a = first(phase_ix)
        rhows = reference_densities(sys)[a]
        rock_density = state.PolymerRockDensity[cell]
        dps = state.DeadPoreSpace[cell]
        bulk_vol = state.BulkVolume[cell]
        rock_vol = bulk_vol - pore_vol
        ads = state.AdsorbedPolymerConcentration[cell]
        mass = v*(1.0 - dps) + rhows*rock_vol*rock_density*ads
    end
    return mass
end

function Jutul.get_dependencies(tracer::PolymerTracer, model)
    if model_or_domain_is_well(model)
        out = Symbol[]
    else
        out = [
            :AdsorbedPolymerConcentration,
            :PolymerConcentration,
            :BulkVolume,
            :DeadPoreSpace,
            :PolymerRockDensity
        ]
    end
    return out
end

function tracer_total_mass_flux(tracer::PolymerTracer, model, state, phase_mass_fluxes, index, upw, C = state.TracerConcentrations, T = eltype(C))
    phase = tracer.ix
    q_ph = phase_mass_fluxes[phase]
    M = state.EffectivePolymerViscosityMultipliers
    F = cell -> C[index, cell]/M[2, cell]
    C_iface = upwind(upw, F, q_ph)
    return C_iface*q_ph
end

struct MaxPolymerConcentration <: Jutul.ScalarVariable end

function Jutul.update_parameter_before_step!(s_max, ::MaxPolymerConcentration, storage, model, dt, forces)
    s = storage.state.PolymerConcentration
    JutulDarcy.update_max_hysteresis_value!(s_max, s)
    return s_max
end

struct DeadPoreSpace{R} <: Jutul.ScalarVariable
    default_dps::Vector{R}
end

function Jutul.default_values(model, var::DeadPoreSpace)
    ne = number_of_entities(model, var)
    dps = var.default_dps
    @assert length(dps) == ne "Number of default values does not match number of entities"
    return copy(dps)
end

struct PolymerRockDensity{R} <: Jutul.ScalarVariable
    default_rho::Vector{R}
end

function Jutul.default_values(model, var::PolymerRockDensity)
    # Note: In principle this could be different than the regular rock density
    # used in thermal. We duplicate it just in case.
    ne = number_of_entities(model, var)
    dps = var.default_rho
    @assert length(dps) == ne "Number of default values does not match number of entities"
    return copy(dps)
end

struct PolymerConcentration <: Jutul.ScalarVariable
    tracer_ix::Int
end

Jutul.@jutul_secondary function update_polymer_concentration!(vals, def::PolymerConcentration, model, TracerConcentrations, Saturations, ix)
    a = first(JutulDarcy.phase_indices(model.system))
    for i in ix
        sw = Saturations[a, i]
        if sw > 0
            v = TracerConcentrations[def.tracer_ix, i]
        else
            v = 0.0
        end
        vals[i] = v
    end
    return vals
end

struct AdsorbedPolymerConcentration{R, A, P} <: Jutul.ScalarVariable
    plyrock::R
    plyads::A
    regions::P
    function AdsorbedPolymerConcentration(plyrock, plyads, satnum = nothing)
        @assert length(plyrock) == length(plyads)
        JutulDarcy.check_regions(satnum, length(plyrock))
        return new{typeof(plyrock), typeof(plyads), typeof(satnum)}(plyrock, plyads, satnum)
    end
end

Jutul.@jutul_secondary function update_polymer_adsorption!(vals, def::AdsorbedPolymerConcentration, model, PolymerConcentration, MaxPolymerConcentration, ix)
    for cell in ix
        reg = region(def.regions, cell)
        plyrock = table_by_region(def.plyrock, reg)
        plyads = table_by_region(def.plyads, reg)
        cp = PolymerConcentration[cell]
        if plyrock.desorption
            eval_c = cp
        else
            cpmax = MaxPolymerConcentration[cell]
            eval_c = max(cp, cpmax)
        end
        vals[cell] = plyads(eval_c)
    end
    return vals
end

struct FullyMixedPolymerViscosityMultiplier{T, C, R} <: Jutul.JutulVariables
    mixed_polymer_viscosity::T
    max_mixed_multiplier::C
    regions::R
    function FullyMixedPolymerViscosityMultiplier(mixed_polymer_viscosity::T, max_concentration, viscosity_regions::R = nothing) where {T, R}
        JutulDarcy.check_regions(viscosity_regions, length(mixed_polymer_viscosity))
        max_mixed_multiplier = map((f, c) -> f(c), mixed_polymer_viscosity, max_concentration)
        return new{T, typeof(max_mixed_multiplier), R}(mixed_polymer_viscosity, max_mixed_multiplier, viscosity_regions)
    end
end

Jutul.degrees_of_freedom_per_entity(model, ::FullyMixedPolymerViscosityMultiplier) = 2


Jutul.@jutul_secondary function update_polymer_viscosity!(vals, def::FullyMixedPolymerViscosityMultiplier, model, PolymerConcentration, ix)
    for cell in ix
        c = PolymerConcentration[cell]
        r = region(def.regions, cell)
        F = table_by_region(def.mixed_polymer_viscosity, r)
        vals[1, cell] = F(c)
        vals[2, cell] = def.max_mixed_multiplier[r]
    end
    return vals
end

struct EffectivePolymerViscosityMultipliers{F, R} <: Jutul.VectorVariables
    mixpar::F
    max_concentration::F
    rrf::Vector{F}
    ads_max::Vector{F}
    regions::R
    function EffectivePolymerViscosityMultipliers(;
            max_concentration::F,
            mixpar::F,
            rrf::Vector{F},
            ads_max::Vector{F},
            regions = nothing,
        ) where F
        JutulDarcy.check_regions(regions, length(rrf))

        return new{F, typeof(regions)}(
            mixpar,
            max_concentration,
            rrf,
            ads_max,
            regions,
        )
    end
end

Jutul.degrees_of_freedom_per_entity(model, ::EffectivePolymerViscosityMultipliers) = 2

Jutul.@jutul_secondary function update_polymer_multipliers!(vals, def::EffectivePolymerViscosityMultipliers, model, FullyMixedPolymerViscosityMultiplier, PolymerConcentration, AdsorbedPolymerConcentration, ix)
    water_ind = first(JutulDarcy.phase_indices(model.system))
    for cell in ix
        c = PolymerConcentration[cell]
        ads = AdsorbedPolymerConcentration[cell]
        mult = FullyMixedPolymerViscosityMultiplier[1, cell]
        mult_max = FullyMixedPolymerViscosityMultiplier[2, cell]
        mu_w_mult, mu_p_mult = polymer_multipliers(def, c, mult, mult_max, ads, cell)
        vals[1, cell] = mu_w_mult
        vals[2, cell] = mu_p_mult
    end
    return vals
end

function polymer_multipliers(def::EffectivePolymerViscosityMultipliers, c, mult, mult_max, ads, cell)
    mixpar = def.mixpar
    c_max = def.max_concentration
    reg_rock = region(def.regions, cell)
    rrf = table_by_region(def.rrf, reg_rock)
    ads_max = table_by_region(def.ads_max, reg_rock)
    ads_adjustment = 1.0 + (rrf - 1.0)*ads/ads_max
    c_norm = c/c_max

    α = mult_max^(1.0 - mixpar)
    β = 1.0/(1.0 - c_norm + c_norm/α)
    mu_w_mult = ads_adjustment*β*mult^mixpar
    mu_p_mult = α + (1.0-α)*c_norm
    return (mu_w_mult, mu_p_mult)
end

struct PolymerAdjustedViscosities <: JutulDarcy.PhaseVariables
end

Jutul.@jutul_secondary function update_mixed_polymer_viscosity!(vals, def::PolymerAdjustedViscosities, model, BasePhaseViscosities, EffectivePolymerViscosityMultipliers, ix)
    water_ind = first(JutulDarcy.phase_indices(model.system))
    nph = size(vals, 1)

    for cell in ix
        for ph in 1:nph
            val_phase = BasePhaseViscosities[ph, cell]
            if ph == water_ind
                val_phase *= EffectivePolymerViscosityMultipliers[1, cell]
            end
            vals[ph, cell] = val_phase
        end
    end
    return vals
end

function set_polymer_model!(outer_model::MultiModel, datafile; is_well = false)
    rmodel = reservoir_model(outer_model)
    set_polymer_model!(rmodel, datafile, is_well = false)
    for (mname, model) in pairs(outer_model.models)
        if JutulDarcy.model_or_domain_is_well(model)
            set_polymer_model!(model, datafile, is_well = true)
        end
    end
    return outer_model
end

function set_polymer_model!(model::SimulationModel, datafile; is_well = false)
    haskey(datafile["RUNSPEC"], "POLYMER") || throw(ArgumentError("POLYMER keyword not found in RUNSPEC section of datafile"))
    plyvisc = map(plyvisc_table, datafile["PROPS"]["PLYVISC"])
    plyrock = map(plyrock_table, datafile["PROPS"]["PLYROCK"])
    plyads = map(plyads_table, datafile["PROPS"]["PLYADS"])
    plmixpar = plmixpar_parse(datafile["PROPS"]["PLMIXPAR"])
    plymax = plymax_table(datafile["PROPS"]["PLYMAX"])

    if is_well
        satnum = ones(Int, number_of_cells(model.domain))
    else
        # Do the checks only for reservoir model
        reservoir = reservoir_domain(model)
        pvtnum = reservoir[:pvtnum]
        satnum = reservoir[:satnum]
        satnum_max = maximum(satnum)
        @assert maximum(pvtnum) <= length(plyvisc) "PVTNUM values exceed number of PLYVISC tables"
        @assert satnum_max <= length(plyrock) "SATNUM values exceed number of PLYROCK tables"
        @assert satnum_max <= length(plyads) "SATNUM values exceed number of PLYADS tables"
    end

    # Handle polymer concentration
    tracer_ix = findfirst(x -> isa(x, PolymerTracer), model.equations[:tracers].flux_type.tracers)
    @assert !isnothing(tracer_ix) "No polymer tracer found in model"
    set_secondary_variables!(model,
        PolymerConcentration = PolymerConcentration(tracer_ix),
        FullyMixedPolymerViscosityMultiplier = FullyMixedPolymerViscosityMultiplier(plyvisc, plymax, satnum),
    )
    if !is_well
        set_secondary_variables!(model,
            AdsorbedPolymerConcentration = AdsorbedPolymerConcentration(plyrock, plyads, satnum)
        )
        set_parameters!(model,
            MaxPolymerConcentration = MaxPolymerConcentration(),
            DeadPoreSpace = DeadPoreSpace(map(x -> x.dps, plyrock[satnum])),
            PolymerRockDensity = PolymerRockDensity(map(x -> x.rho_rock, plyrock[satnum]))
        )
        vars = Jutul.get_variables_by_type(model, :all)
        if !haskey(vars, :BulkVolume)
            set_parameters!(model,
                BulkVolume = JutulDarcy.BulkVolume(),
            )
        end

        prms = Jutul.get_variables_by_type(model, :parameters)
        svars = Jutul.get_variables_by_type(model, :secondary)
        if haskey(prms, :PhaseViscosities)
            mu = prms[:PhaseViscosities]
            prms[:BasePhaseViscosities] = mu
            delete!(prms, :PhaseViscosities)
        else
            @assert haskey(svars, :PhaseViscosities) "PhaseViscosities not found in secondary variables or parameters, cannot setup polymer model."
            mu = svars[:PhaseViscosities]
            svars[:BasePhaseViscosities] = mu
            delete!(svars, :PhaseViscosities)
        end

        svars[:EffectivePolymerViscosityMultipliers] = EffectivePolymerViscosityMultipliers(
            mixpar = plmixpar,
            rrf = map(x -> x.rrf, plyrock),
            ads_max = map(x -> x.ads_max, plyrock),
            max_concentration = plymax[1],
            regions = pvtnum
        )
        svars[:PhaseViscosities] = PolymerAdjustedViscosities()
    end
    return model
end

function plyvisc_table(tab)
    @assert tab[1, 1] ≈ 0.0 "First entry in column 1 of PLYVISC table must be zero"
    @assert tab[1, 2] ≈ 1.0 "First entry in column 2 of PLYVISC table must be one"
    return get_1d_interpolator(tab[:, 1], tab[:, 2])
end

function plyrock_table(tab)
    dps, rrf, rho_rock, ix, ads_max = tab
    @assert dps >= 0.0 "Dead pore space must be non-negative"
    @assert dps <= 1.0 "Dead pore space must be less than or equal to one"
    @assert rrf >= 1.0 "Residual resistance factor must be greater than or equal to one"
    @assert rho_rock > 0.0 "Rock density must be positive"
    @assert ix in (1, 2) "Adsorption type must be one (desorption) or two (no desorption)"
    @assert ads_max >= 0.0 "Adsorption maximum must be non-negative"

    return (dps = dps, rrf = rrf, rho_rock = rho_rock, desorption = ix == 1, ads_max = ads_max)
end

function plyads_table(tab)
    return get_1d_interpolator(tab[:, 1], tab[:, 2])
end

function plmixpar_parse(tab)
    if length(tab) > 1 
        jutul_message("PLYMIXPAR", "Table must have exactly one row (NPLMIX > 1 not supported). Taking first value.", color = :yellow)
    end
    return only(tab[1])
end

function plymax_table(tab)
    if length(tab) > 1 
        jutul_message("PLYMAX", "Table must have exactly one row (NPLMIX > 1 not supported). Taking first value.", color = :yellow)
    end
    v1, v2 = tab[1]
    return return (v1, v2)
end

function tracer_scale(model, tracer::PolymerTracer)
    # Water density
    return 1000.0
end
