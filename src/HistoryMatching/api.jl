##############################
#             API            #
##############################

function history_match_objective(case::JutulCase, arg...; is_global::Bool = true, kwarg...)
    hm = HistoryMatch(case, arg...; kwarg...)
    if is_global
        return GlobalHistoryMatchObjective(hm)
    else
        ismissing(hm.periods) || error("SumHistoryMatchObjective does not support periods.")
        return SumHistoryMatchObjective(hm)
    end
end

function history_match_objective(hm::HistoryMatch; is_global::Bool = true)
    if is_global
        return GlobalHistoryMatchObjective(hm)
    else
        return SumHistoryMatchObjective(hm)
    end
end

function evaluate_match(obj::HistoryMatchObjective, result::ReservoirSimResult)
    hm = obj.match
    c = hm.case
    model = c.model
    states = result.states
    timesteps = c.dt
    all_forces = c.forces
    return Jutul.evaluate_objective(obj, model, result.result.states, timesteps, all_forces)
end

function match_well!(hm::HistoryMatch, name::Union{String, Symbol}, quantity::Union{String, Symbol};
        weight::Union{Float64, Vector{Float64}} = 1.0,
        is_injector::Bool,
        data = missing,
        t = missing,
        scale = missing,
        is_sum_objective::Bool = false
    )
    is_cumulative = false
    rmodel = reservoir_model(hm.case)
    sys = rmodel.system
    rhows = phase_reference_density(sys, JutulDarcy.AqueousPhase())
    rhoos = phase_reference_density(sys, JutulDarcy.LiquidPhase())
    rhogs = phase_reference_density(sys, JutulDarcy.VaporPhase())
    rhols = 0.5*(rhoos + rhows)

    day = si_unit(:day)
    # Rates
    grat_scale = day/rhogs
    orat_scale = day/rhoos
    wrat_scale = day/rhows
    lrat_scale = day/rhols
    bhp_scale = 1.0/(100*si_unit(:bar))
    # Cumulative production
    gtotal_scale = 1.0./rhogs
    ototal_scale = 1.0./rhoos
    wtotal_scale = 1.0./rhows
    ltotal_scale = 1.0./rhols

    if is_injector
        if quantity isa Symbol
            if quantity == :rate
                quantity = "RATE"
            elseif quantity == :grat
                quantity = "WGIR"
            elseif quantity == :orat
                quantity = "WOIR"
            elseif quantity == :wrat
                quantity = "WWIR"
            elseif quantity == :bhp
                quantity = "WBHP"
            else
                error("Unsupported quantity '$quantity' for injector well match.")
            end
        end
        quantity = uppercase(quantity)
        if quantity == "RATE"
            dest = hm.injector_rate
            dscale = lrat_scale
        elseif quantity == "WGIR"
            dest = hm.injector_grat
            dscale = grat_scale
        elseif quantity == "WOIR"
            dest = hm.injector_orat
            dscale = orat_scale
        elseif quantity == "WWIR"
            dest = hm.injector_wrat
            dscale = wrat_scale
        elseif quantity == "WBHP"
            dest = hm.injector_bhp
            dscale = bhp_scale
        else
            error("Unsupported quantity '$quantity' for injector well match.")
        end
    else
        if quantity isa Symbol
            if quantity == :rate
                quantity = "RATE"
            elseif quantity == :bhp
                quantity = "WBHP"
            elseif quantity == :wrat
                quantity = "WWPR"
            elseif quantity == :grat
                quantity = "WGPR"
            elseif quantity == :orat
                quantity = "WOPR"
            elseif quantity == :lrat
                quantity = "WLPR"
            elseif quantity == :cumulative_oil
                quantity = "WOPT"
            elseif quantity == :cumulative_gas
                quantity = "WGPT"
            elseif quantity == :cumulative_water
                quantity = "WWPT"
            elseif quantity == :cumulative_liquid
                quantity = "WLPT"
            else
                error("Unsupported quantity '$quantity' for producer well match.")
            end
        end
        quantity = uppercase(quantity)
        if quantity == "RATE"
            dest = hm.producer_rate
            dscale = lrat_scale
        elseif quantity == "WBHP"
            dest = hm.producer_bhp
            dscale = bhp_scale
        elseif quantity == "WGPR"
            dest = hm.producer_grat
            dscale = grat_scale
        elseif quantity == "WOPR"
            dest = hm.producer_orat
            dscale = orat_scale
        elseif quantity == "WLPR"
            dest = hm.producer_lrat
            dscale = lrat_scale
        elseif quantity == "WWPR"
            dest = hm.producer_wrat
            dscale = wrat_scale
        elseif quantity == "WOPT"
            dest = hm.producer_cumulative_oil
            is_cumulative = true
            dscale = gtotal_scale
        elseif quantity == "WGPT"
            dest = hm.producer_cumulative_gas
            is_cumulative = true
            dscale = ototal_scale
        elseif quantity == "WWPT"
            dest = hm.producer_cumulative_water
            is_cumulative = true
            dscale = wtotal_scale
        elseif quantity == "WLPT"
            dest = hm.producer_cumulative_liquid
            is_cumulative = true
            dscale = ltotal_scale
        else
            error("Unsupported quantity '$quantity' for producer well match.")
        end
    end
    if is_sum_objective && is_cumulative
        error("Cumulative well matches are not supported for SumHistoryMatchObjective.")
    end
    data = get_well_data(hm, name, quantity, data, t)
    if weight isa Vector
        length(weight) == length(case.dt) || error("Weight vector length must match number of simulation report steps in case ($(length(case.dt))).")
    end
    if ismissing(scale)
        scale = dscale
    end
    wm = WellMatch(Symbol(name), quantity, weight, scale, is_injector, data)
    # total rates can be computed from other rates...
    push!(dest, wm)
    return hm
end

function phase_reference_density(sys, phase)
    if phase in JutulDarcy.get_phases(sys)
        idx = JutulDarcy.phase_index(sys, VaporPhase())
        rhos = JutulDarcy.reference_densities(sys)
        val = rhos[idx]
    else
        val = 1000.0
    end
    return val
end

function match_injectors!(hm::HistoryMatch, quantity::Union{String, Symbol}, wells = get_injectors(hm); kwarg...)
    if wells isa Symbol || wells isa String
        wells = [wells]
    end
    for w in wells
        match_well!(hm, w, quantity; is_injector = true, kwarg...)
    end
    return hm
end

function match_producers!(hm::HistoryMatch, quantity::Union{String, Symbol}, wells = get_producers(hm); kwarg...)
    if wells isa Symbol || wells isa String
        wells = [wells]
    end
    for w in wells
        match_well!(hm, w, quantity; is_injector = false, kwarg...)
    end
    return hm
end

match_injectors!(obj::HistoryMatchObjective, arg...; kwarg...) = match_injectors!(obj.match, arg...; kwarg...)
match_producers!(obj::HistoryMatchObjective, arg...; kwarg...) = match_producers!(obj.match, arg...; kwarg...)
match_well!(obj::GlobalHistoryMatchObjective, arg...; kwarg...) = match_well!(obj.match, arg...; kwarg...)
match_well!(obj::SumHistoryMatchObjective, arg...; kwarg...) = match_well!(obj.match, arg...; kwarg..., is_sum_objective = true)
