import JutulDarcy: region, table_by_region

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

function tracer_total_mass_outer(tracer::PolymerTracer, model, state, concentration, resident_mass_density, vol, cell, index)
    if resident_mass_density <= TRACER_TOL
        v = vol*Jutul.replace_value(concentration, 0.0)
    else
        v = tracer_total_mass(tracer, model, state, concentration, resident_mass_density, vol, cell, index)
    end
    if JutulDarcy.model_or_domain_is_well(model)
        mass = v
    else
        rock_density = state.PolymerRockDensity[cell]
        dps = state.DeadPoreSpace[cell]
        bulk_vol = state.BulkVolume[cell]
        ads = state.AdsorbedPolymerConcentration[cell]
        mass = v*(1 - dps) + (vol - bulk_vol)*rock_density*ads
    end
    return mass
end

function Jutul.get_dependencies(tracer::PolymerTracer, model)
    if JutulDarcy.model_or_domain_is_well(model)
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
    M = state.PolymerViscosityMultipliers
    F = cell -> C[index, cell]/M[2, cell]
    C_iface = JutulDarcy.upwind(upw, F, q_ph)
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

struct PolymerViscosityMultipliers{T, R, RREG, VREG} <: Jutul.VectorVariables
    mixed_polymer_viscosity::T
    mixpar::R
    max_concentration::R
    rrf::Vector{R}
    ads_max::Vector{R}
    rock_regions::RREG
    viscosity_regions::VREG
    function PolymerViscosityMultipliers(
        mixed_polymer_viscosity::T;
        max_concentration::R,
        mixpar::R,
        rrf::Vector{R},
        ads_max::Vector{R},
        rock_regions = nothing,
        viscosity_regions = nothing
        ) where {T, R}
        JutulDarcy.check_regions(viscosity_regions, length(mixed_polymer_viscosity))
        JutulDarcy.check_regions(rock_regions, length(rrf))

        return new{T, R, typeof(rock_regions), typeof(viscosity_regions)}(
            mixed_polymer_viscosity,
            mixpar,
            max_concentration,
            rrf,
            ads_max,
            rock_regions,
            viscosity_regions
        )
    end
end

Jutul.degrees_of_freedom_per_entity(model, ::PolymerViscosityMultipliers) = 3

Jutul.@jutul_secondary function update_polymer_multipliers!(vals, def::PolymerViscosityMultipliers, model, PolymerConcentration, AdsorbedPolymerConcentration, ix)
    water_ind = first(JutulDarcy.phase_indices(model.system))
    for cell in ix
        c = PolymerConcentration[cell]
        ads = AdsorbedPolymerConcentration[cell]

        mu_w_mult, mu_p_mult, cmult = polymer_multipliers(def, c, ads, cell)
        vals[1, cell] = mu_w_mult
        vals[2, cell] = mu_p_mult
        vals[3, cell] = cmult
    end
    return vals
end

function polymer_multipliers(def::PolymerViscosityMultipliers, c, ads, cell)
    mixpar = def.mixpar
    c_max = def.max_concentration
    reg_rock = region(def.rock_regions, cell)
    rrf = table_by_region(def.rrf, reg_rock)
    ads_max = table_by_region(def.ads_max, reg_rock)
    reg_visc = region(def.viscosity_regions, cell)
    mu = table_by_region(def.mixed_polymer_viscosity, reg_visc)
    ads_adjustment = 1.0 + (rrf - 1.0)*ads/ads_max
    c_norm = c/c_max
    mult = mu(c)
    mult_max = mu(c_max)

    α = mult_max^(1.0 - mixpar)
    β = 1.0/(1.0 - c_norm + c_norm/α)
    mu_w_mult = ads_adjustment*β*mult^mixpar
    mu_p_mult = α + (1.0-α)*c_norm
    return (mu_w_mult, mu_p_mult, mult)
end

struct PolymerAdjustedViscosities <: JutulDarcy.PhaseVariables
end

Jutul.@jutul_secondary function update_mixed_polymer_viscosity!(vals, def::PolymerAdjustedViscosities, model, BasePhaseViscosities, PolymerViscosityMultipliers, ix)
    water_ind = first(JutulDarcy.phase_indices(model.system))
    nph = size(vals, 1)

    for cell in ix
        for ph in 1:nph
            val_phase = BasePhaseViscosities[ph, cell]
            if ph == water_ind
                val_phase *= PolymerViscosityMultipliers[1, cell]
            end
            vals[ph, cell] = val_phase
        end
    end
    return vals
end

function set_polymer_model!(outer_model, datafile)
    Jutul.jutul_message("POLYMER", "Polymer model is in early development. Use with caution.", color = :yellow)
    model = reservoir_model(outer_model)
    reservoir = reservoir_domain(model)
    haskey(datafile["RUNSPEC"], "POLYMER") || throw(ArgumentError("POLYMER keyword not found in RUNSPEC section of datafile"))

    pvtnum = reservoir[:pvtnum]
    satnum = reservoir[:satnum]
    satnum_max = maximum(satnum)

    plyvisc = map(plyvisc_table, datafile["PROPS"]["PLYVISC"])
    @assert maximum(pvtnum) <= length(plyvisc) "PVTNUM values exceed number of PLYVISC tables"
    plyrock = map(plyrock_table, datafile["PROPS"]["PLYROCK"])
    @assert satnum_max <= length(plyrock) "SATNUM values exceed number of PLYROCK tables"
    plyads = map(plyads_table, datafile["PROPS"]["PLYADS"])
    @assert satnum_max <= length(plyads) "SATNUM values exceed number of PLYADS tables"
    plmixpar = plmixpar_parse(datafile["PROPS"]["PLMIXPAR"])
    plymax = plymax_table(datafile["PROPS"]["PLYMAX"])
    # Handle polymer concentration
    tracer_ix = findfirst(x -> isa(x, PolymerTracer), model.equations[:tracers].flux_type.tracers)
    @assert !isnothing(tracer_ix) "No polymer tracer found in model"
    set_secondary_variables!(model,
        PolymerConcentration = PolymerConcentration(tracer_ix),
        # PolymerTotalMass = PolymerTotalMass(plyrock, plyads, satnum),
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

    svars[:PolymerViscosityMultipliers] = PolymerViscosityMultipliers(
        plyvisc,
        mixpar = plmixpar,
        rrf = map(x -> x.rrf, plyrock),
        ads_max = map(x -> x.ads_max, plyrock),
        max_concentration = plymax[1],
        rock_regions = pvtnum,
        viscosity_regions = satnum
    )
    svars[:PhaseViscosities] = PolymerAdjustedViscosities()
    return outer_model
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
