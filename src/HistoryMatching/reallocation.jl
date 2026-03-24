function reallocate_oil_production_rates(smry, group;
        start::Int = 1,
        stop::Int = length(smry["TIME"].seconds),
        oil_error = 0.075,
        weights = rand(length(group))
    )
    reallocation = weights ./ sum(weights)

    smry = deepcopy(smry)
    wdata = smry["VALUES"]["WELLS"]

    function redistribute_linear(oil_rates, liquid_rates, oil_to_distribute, weights)
        return min.(oil_rates + weights.*oil_to_distribute, liquid_rates)
    end


    for step in start:stop
        orates = map(w -> wdata[w]["WOPR"][step], group)
        lrates = map(w -> wdata[w]["WLPR"][step], group)

        total_oil_rate = sum(orates)
        new_orates = orates.*(1.0 - oil_error)
        excess_oil = oil_error*total_oil_rate
        @assert sum(new_orates) + excess_oil ≈ sum(orates)

        remaining_oil = excess_oil
        it = 1
        w = copy(reallocation)
        while remaining_oil > 1e-12
            w[new_orates .≈ lrates] .= 0.0
            w = w./sum(w)
            new_orates = redistribute_linear(new_orates, lrates, remaining_oil, w)

            remaining_oil = total_oil_rate - sum(new_orates)
            it += 1
        end
        for (i, w) in enumerate(group)
            new_oil_rate = new_orates[i]
            new_water_rate = wdata[w]["WLPR"][step] - new_oil_rate
            wdata[w]["WOPR"][step] = new_oil_rate
            wdata[w]["WWPR"][step] = new_water_rate
            # Now, update the cumulative oil and water rates
            if step == 1
                prev_w = 0.0
                prev_o = 0.0
                t_prev = 0.0
            else
                t_prev = smry["TIME"].seconds[step-1]
                prev_w = wdata[w]["WWPT"][step-1]
                prev_o = wdata[w]["WOPT"][step-1]
            end
            dt = smry["TIME"].seconds[step] - t_prev
            wdata[w]["WWPT"][step] = prev_w + wdata[w]["WWPR"][step]*dt
            wdata[w]["WOPT"][step] = prev_o + wdata[w]["WOPR"][step]*dt
        end
    end
    return smry
end

