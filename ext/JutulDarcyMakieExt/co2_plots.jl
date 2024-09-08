function JutulDarcy.plot_co2_inventory(t, inventory; plot_type = :stack)
    getd(x) = map(i -> i[x], inventory)
    function get_label(k)
        if k == :dissolved
            label = "Dissolved"
        elseif k == :mobile
            label = "Mobile"
        elseif k == :residual
            label = "Residual"
        elseif k == :inside_total
            label = "Inside region"
        elseif k == :outside_total
            label = "Outside region"
        elseif k == :outside_domain
            label = "Outside domain"
        else
            label = "$k"
        end
    end
    is_subset = any(x -> x[:outside_total] > 1e-3, inventory)
    if is_subset
        targets = [:dissolved, :residual, :mobile, :outside_total, :outside_domain]
    else
        targets = [:dissolved, :residual, :mobile, :outside_domain]
    end
    if eltype(t) == Float64
        x = t./3.1556952e7
        xl = "Years"
    else
        x = t
        xl = "Step"
    end
    n = length(inventory)
    fig = Figure()
    ax = Axis(fig[1, 1], ylabel = "CO2 mass (kg)", xlabel = xl)
    cake_data = Float64[]
    if plot_type == :stack
        prev = zeros(n)
        for k in targets
            val = getd(k)
            push!(cake_data, val[end])
            lower = prev
            upper = lower .+ val
            band!(ax, x, lower, upper, label = get_label(k))
            @. prev = upper
        end
    else
        @assert plot_type == :lines
        for k in targets
            val = getd(k)
            push!(cake_data, val[end])
            lines!(ax, x, lower, upper, label = get_label(k))
        end
    end
    axislegend(position = :lt)
    return fig
end

function JutulDarcy.plot_co2_inventory(inventory::Vector)
    JutulDarcy.plot_co2_inventory(eachindex(inventory), inventory)
end
