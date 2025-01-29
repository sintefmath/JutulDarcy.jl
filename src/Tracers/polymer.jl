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

tracer_phase_indices(t::PolymerTracer) = (t.phase_index, )

function tracer_total_mass_outer(tracer::PolymerTracer, model, state, concentration, resident_mass_density, vol, cell, index)
    error()
    water_mass = resident_mass
    # flowing = water_mass*
    # return concentration*resident_mass
end



struct PolymerConcentration <: Jutul.ScalarVariable
    tracer_ix::Int
end

struct MaxPolymerConcentration <: Jutul.ScalarVariable end

function Jutul.update_parameter_before_step!(s_max, ::MaxPolymerConcentration, storage, model, dt, forces)
    s = storage.state.PolymerConcentration
    update_max_hysteresis_value!(s_max, s)
    return s_max
end

Jutul.@jutul_secondary function update_polymer_concentration(vals, def::PolymerConcentration, model, TracerConcentrations, Saturations, ix)
    a = first(phase_indices(model.system))
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


struct PolymerTotalMass{R, A, V} <: Jutul.ScalarVariable
    plyrock::R
    plyads::A
    regions::V
    function PolymerTotalMass(plyrock, plyads, satnum = nothing)
        @assert length(plyrock) == length(plyads)
        JutulDarcy.check_regions(satnum, length(plyrock))
        return new{typeof(plyrock), typeof(plyads), typeof(satnum)}(plyrock, plyads, satnum)
    end
end

function set_polymer_model!(model, datafile)
    model = reservoir_model(model)
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
        PolymerTotalMass = PolymerTotalMass(plyrock, plyads, satnum),
    )
    set_parameters!(model,
        MaxPolymerConcentration = MaxPolymerConcentration(),
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
        prms[:BasePhaseViscosities] = plyvisc
        mu = prms[:PhaseViscosities]
        delete!(prms, :PhaseViscosities)
    else
        @assert haskey(svars, :PhaseViscosities) "PhaseViscosities not found in secondary variables or parameters, cannot setup polymer model."
        svars[:BasePhaseViscosities] = plyvisc
        mu = svars[:PhaseViscosities]
        delete!(svars, :PhaseViscosities)
    end

    # Set up to manage the accumulation term

    # Set up to manage the adsorption term

    # Set up to manage the viscosity modification


    error()
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
