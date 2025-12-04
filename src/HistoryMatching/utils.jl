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
    if mint > 0.0
        jutul_message("HistoryMatch", "Minimum time in provided data for well $name, $quantity is greater than zero ($(Jutul.get_tstr(mint))). Constant extrapolation to time zero will be performed.")
    end
    if maxt > sum(hm.case.dt)
        jutul_message("HistoryMatch", "Maximum time in provided data for well $name, $quantity exceeds total simulation time ($(Jutul.get_tstr(sum(hm.case.dt)))). Time must be given in seconds.")
    end
    minval = minimum(response)
    if minval <= 0.0
        jutul_message("HistoryMatch", "Minimum value in provided data for well $name, $quantity is negative ($minval). All well responses are assumed to be non-negative and into/out of reservoir is set via injector/producer designation.")
    end
    return Jutul.get_1d_interpolator(time, response, constant_dx = false, static = false)
end
