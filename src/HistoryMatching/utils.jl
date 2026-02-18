function set_report_step_and_weights_for_period!(period::MatchPeriod, dt::Vector{Float64}, baseweight)
    start, stop = period.start, period.stop
    idxs = period.step_idx
    fracs = period.weights
    length(idxs) == 0 || error("Period step indices have already been set.")
    length(fracs) == 0 || error("Period step weights have already been set.")
    t_start = 0.0
    for (i, dt_i) in enumerate(dt)
        t_stop = t_start + dt_i
        if start <= t_stop && t_start <= stop
            w_i = get_period_weight(period, baseweight, t_start, t_stop)
            push!(idxs, i)
        else
            w_i = 0.0
        end
        push!(fracs, w_i)
        t_start = t_stop
    end
    return period
end


function get_period_weight(p::MatchPeriod, baseweight, start, stop)
    Δstep = stop - start
    Δoverlap = min(p.stop, stop) - max(p.start, start)
    return baseweight*(Δoverlap/Δstep)
end

function setup_periods(case, periods, period_weights, normalized_periods)
    if eltype(periods) == MatchPeriod
        ismissing(period_weights) || error("Cannot specify period_weights when periods are MatchPeriod.")
        ismissing(normalized_periods) || error("Cannot specify normalized_periods when periods are MatchPeriod.")
    else
        if !ismissing(normalized_periods)
            ismissing(period_weights) || error("Cannot specify period_weights when normalized_periods are given.")
            total_time = sum(case.dt)
            periods = Float64[]
            for t in normalized_periods
                push!(periods, t*total_time)
            end
        end
    end
    if eltype(periods) == MatchPeriod
        # Do nothing
    elseif !ismissing(periods)
        new_periods = MatchPeriod[]
        if eltype(periods)<:Real
            prev = 0.0
            for t in periods
                # w = period_weights[i]
                push!(new_periods, MatchPeriod(prev, t))
                prev = t
            end
            periods = new_periods
        else
            for tup in periods
                length(tup) == 2 || error("Each period tuple must have exactly two elements (start, stop).")
                prev, next = tup
                # w = period_weights[i]
                push!(new_periods, MatchPeriod(prev, next))
            end
            periods = new_periods
        end
    end
    if !ismissing(periods)
        if ismissing(period_weights)
            period_weights = ones(length(periods))
        end
        for (i, p) in enumerate(periods)
            p.start < p.stop || error("In period $i, start time must be less than stop time.")
            if i > 1
                periods[i-1].stop <= p.start || error("Periods must be non-overlapping and in increasing order.")
            end
            set_report_step_and_weights_for_period!(p, case.dt, period_weights[i])
            @assert length(p.step_idx) > 0
            @assert length(p.weights) == length(case.dt)
        end
    end
    return periods
end

function get_injectors(hm::HistoryMatch)
    return get_wells_by_control_type(hm, :injector)
end

function get_producers(hm::HistoryMatch)
    return get_wells_by_control_type(hm, :producer)
end

function get_wells(hm::HistoryMatch)
    return get_wells(hm.case)
end

function get_wells(c::JutulCase)
    return collect(keys(JutulDarcy.get_model_wells(c.model)))
end

function get_wells_by_control_type(hm::HistoryMatch, control_type::Symbol)
    if control_type == :injector
        ctype = InjectorControl
    elseif control_type == :producer
        ctype = ProducerControl
    elseif control_type == :disabled
        ctype = DisabledControl
    else
        error("Unknown control type '$control_type'.")
    end
    wells = get_wells(hm)
    c = hm.case
    forces = c.forces
    if !(forces isa AbstractVector)
        forces = [forces]
    end
    output_wells = Symbol[]
    for f in forces
        if length(wells) == 0
            break
        end
        for w in wells
            c = f[:Facility].control[w]
            if c isa ctype
                wells = setdiff(wells, [w])
                push!(output_wells, w)
            end
        end
    end
    return unique!(output_wells)
end

function flexible_getindex(x, key::Symbol)
    if haskey(x, key)
        return x[key]
    elseif haskey(x, "$key")
        return x["$key"]
    else
        return missing
    end
end

function flexible_getindex(x, key::String)
    if haskey(x, key)
        return x[key]
    elseif haskey(x, Symbol(key))
        return x[Symbol(key)]
    else
        return missing
    end
end

function get_well_data(hm::HistoryMatch, name, quantity, data, t)
    if ismissing(data)
        ismissing(t) || error("t was provided, but not data for well $name, $quantity.")
        !ismissing(hm.summary) || error("No summary data available in HistoryMatch. You must provide data and time vector t to match $quantity.")
        if haskey(hm.summary, "VALUES")
            # Jutul style summary object
            smry_data = hm.summary["VALUES"]
            time = hm.summary["TIME"].seconds
            wdata = flexible_getindex(smry_data["WELLS"], name)
            response = flexible_getindex(wdata, quantity)
            if ismissing(response)
                error("Summary data for well $name, quantity $quantity not found in summary object.")
            end
        else
            error("Unsupported summary data format. Cannot extract well data for $name, $quantity.")
        end
    else
        !ismissing(t) || error("If data is provided, time vector t must also be provided.")
        length(data) == length(t) || error("Data and time vector t must have the same length.")
        all(diff(t) .> 0) || error("Time vector t must be strictly increasing.")
        time = t
        response = data
    end
    mint, maxt = extrema(time)
    if mint > 1.01*hm.case.dt[1]
        jutul_message("HistoryMatch", "Minimum time in provided data for well $name, $quantity is greater than zero ($(Jutul.get_tstr(mint))). Constant extrapolation to time zero will be performed.")
    end
    total_case_t = sum(hm.case.dt)
    if maxt > sum(hm.case.dt)
        jutul_message("HistoryMatch", "Maximum time in provided data for well $name, $quantity ($(Jutul.get_tstr(maxt))) exceeds total simulation time ($(Jutul.get_tstr(total_case_t))). Time must be given in seconds.")
    end
    minval = minimum(response)
    if minval < -1e-12
        jutul_message("HistoryMatch", "Minimum value in provided data for well $name, $quantity is negative ($minval). All well responses are assumed to be non-negative and into/out of reservoir is set via injector/producer designation.")
    end
    return Jutul.get_1d_interpolator(time, response, constant_dx = false, static = false)
end

function mismatch_summary(obj::HistoryMatchObjective, res::ReservoirSimResult, fld::String, threshold = 0.2; kwarg...)
    return mismatch_summary(obj.match.summary, res.summary, fld; kwarg...)
end

function mismatch_summary(summary_ref, summary, fld::String, threshold = 0.2;
        wells = keys(summary["VALUES"]["WELLS"]),
        do_print = 0,
        npts = 1000,
        no_data_threshold = ifelse(fld == "WBHP", si_unit(:atm), 1e-10),
        relative = true,
        type = :mean,
        prefix = ""
    )
    t_ref = summary_ref["TIME"].seconds
    t = summary["TIME"].seconds

    t_eval = range(0.0, stop = t[end], length = npts)
    function resample(t, val)
        return Jutul.get_1d_interpolator(t, val).(t_eval)
    end

    start_ref = summary["TIME"].start_date
    start = summary["TIME"].start_date
    if !(ismissing(start) && ismissing(start_ref))
        if start != start_ref
            jutul_message("HistoryMatch", "Mismatch summary: Start time in reference data ($(start_ref)) does not match start time in simulation summary ($(start)). Assuming equal start dates!")
        end
    end

    mismatch = Dict{String, Float64}()
    bad = String[]
    mismatch_per_step = Dict{String, Vector{Float64}}()
    values_per_step = Dict{String, Vector{Float64}}()
    reference_per_step = Dict{String, Vector{Float64}}()
    for w in wells
        val_ref = abs.(resample(t_ref, summary_ref["VALUES"]["WELLS"][w][fld]))
        val_sim = abs.(resample(t, summary["VALUES"]["WELLS"][w][fld]))

        step_mismatch = zeros(length(val_ref))
        num_ok = 0
        for i in eachindex(val_ref)
            if abs(val_ref[i]) < no_data_threshold || !isfinite(val_ref[i])
                continue
            end
            step_mismatch[i] = abs(val_sim[i] - val_ref[i])
            num_ok += 1
        end
        if type == :mean
            mval = sum(step_mismatch)
            if relative
                mval /= max(sum(abs.(val_ref)), no_data_threshold)
            else
                mval /= max(num_ok, 1)
            end
        elseif type == :max
            if relative
                worst = 0.0
                for i in eachindex(step_mismatch)
                    v = val_sim[i]/val_ref[i]
                    if v > worst
                        worst = v
                    end
                end
                mval = worst
            else
                mval = maximum(step_mismatch)
            end
        else
            error("Unknown mismatch type '$type'. Supported types are :mean and :max.")
        end
        mismatch[w] = mval
        mismatch_per_step[w] = step_mismatch
        values_per_step[w] = val_sim
        reference_per_step[w] = val_ref
        if mval >= threshold
            push!(bad, w)
        end
    end
    prt = Int(do_print)
    if prt > 0
        if prt > 1
            # Potentially very detailed printout
            println("$(prefix)Mismatch summary for field '$fld':")
            if length(bad) > 0
                println("Wells with mismatch above threshold ($threshold): ", join(bad, ", "))
                if prt > 2
                    for w in bad
                        mval = mismatch[w]
                        println("  $w: $mval")
                    end
                end
            else
                println("$(prefix)All wells have mismatch below threshold ($threshold).")
            end
        else
            fmt(x) = round(x, sigdigits=2)
            worstval = 0.0
            worstkey = ""
            bestval = typemax(Float64)
            bestkey = ""
            avg = 0.0
            n = 1
            for (k, v) in pairs(mismatch)
                if v > worstval
                    worstval = v
                    worstkey = k
                end
                if v < bestval && any(x -> x > no_data_threshold, summary["VALUES"]["WELLS"][k][fld])
                    bestval = v
                    bestkey = k
                end
                avg += v
                n += 1
            end
            if relative
                fldk = "rel.$type"
            else
                fldk = type
            end
            print("$prefix$fld $fldk: ")
            if length(bad) > 0
                print("$(length(bad))/$(length(wells)) above $threshold")
            else
                print("All wells below threshold ($threshold).")
            end
            println(" (avg $(fmt(avg/n))), worst: $worstkey: $(fmt(worstval)), best $bestkey: $(fmt(bestval))")
        end
    end
    return Dict(
        "mismatch" => mismatch,
        "bad" => bad,
        "mismatch_per_step" => mismatch_per_step,
        "values_per_step" => values_per_step,
        "reference_per_step" => reference_per_step,
        "seconds" => t_eval,
    )
end
