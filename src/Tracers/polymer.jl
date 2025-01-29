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
    error()
    water_mass = resident_mass
    # TODO: Dependencies
    # flowing = water_mass*
    # return concentration*resident_mass
end


struct MaxPolymerConcentration <: Jutul.ScalarVariable end

function Jutul.update_parameter_before_step!(s_max, ::MaxPolymerConcentration, storage, model, dt, forces)
    s = storage.state.PolymerConcentration
    update_max_hysteresis_value!(s_max, s)
    return s_max
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


struct PolymerAdjustedViscosities{T, R, RREG, VREG} <: JutulDarcy.PhaseVariables
    mixed_polymer_viscosity::T
    mixpar::R
    max_concentration::R
    rrf::Vector{R}
    ads_max::Vector{R}
    rock_regions::RREG
    viscosity_regions::VREG
    function PolymerAdjustedViscosities(
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

Jutul.@jutul_secondary function update_mixed_polymer_viscosity!(vals, def::PolymerAdjustedViscosities, model, PolymerConcentration, BasePhaseViscosities, AdsorbedPolymerConcentration, ix)
    water_ind = first(JutulDarcy.phase_indices(model.system))
    mixpar = def.mixpar
    c_max = def.max_concentration
    nph = size(vals, 1)

    for cell in ix
        reg_rock = region(def.rock_regions, cell)
        rrf = table_by_region(def.rrf, reg_rock)
        ads_max = table_by_region(def.ads_max, reg_rock)
        reg_visc = region(def.viscosity_regions, cell)
        mu = table_by_region(def.mixed_polymer_viscosity, reg_visc)

        c = PolymerConcentration[cell]
        ads = AdsorbedPolymerConcentration[cell]
        ads_adjustment = 1.0 + (rrf - 1.0)*ads/ads_max
        c_norm = c/c_max
        α = mu(c_max)^(1.0 - mixpar)
        β = 1.0/(1.0 - c_norm + c_norm/α)
        mu_w_mult = ads_adjustment*β*mu(c)^mixpar

        for ph in 1:nph
            val_phase = BasePhaseViscosities[ph, cell]
            if ph == water_ind
                val_phase *= mu_w_mult
            end
            vals[ph, cell] = val_phase
        end
    end
    return vals
end

# struct PolymerTotalMass{R, A, V} <: Jutul.ScalarVariable
#     plyrock::R
#     plyads::A
#     regions::V
#     function PolymerTotalMass(plyrock, plyads, satnum = nothing)
#         @assert length(plyrock) == length(plyads)
#         JutulDarcy.check_regions(satnum, length(plyrock))
#         return new{typeof(plyrock), typeof(plyads), typeof(satnum)}(plyrock, plyads, satnum)
#     end
# end

function set_polymer_model!(outer_model, datafile)
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
    @info "???" plyads plyrock plmixpar plymax
    set_secondary_variables!(model,
        PolymerConcentration = PolymerConcentration(tracer_ix),
        # PolymerTotalMass = PolymerTotalMass(plyrock, plyads, satnum),
        AdsorbedPolymerConcentration = AdsorbedPolymerConcentration(plyrock, plyads, satnum)
    )
    set_parameters!(model,
        MaxPolymerConcentration = MaxPolymerConcentration()
    )
    vars = Jutul.get_variables_by_type(model, :all)
    if !haskey(vars, :BulkVolume)
        set_parameters!(model,
            BulkVolume = JutulDarcy.BulkVolume()
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

    svars[:PhaseViscosities] = PolymerAdjustedViscosities(
        plyvisc,
        mixpar = plmixpar,
        rrf = map(x -> x.rrf, plyrock),
        ads_max = map(x -> x.ads_max, plyrock),
        max_concentration = plymax[1],
        rock_regions = pvtnum,
        viscosity_regions = satnum
    )
    # error()
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
