struct CSP11ReportHelper
    specase::Symbol
    points::Tuple{Int, Int}
    regions::Tuple{Vector{Tuple{Int, Float64}}, Vector{Tuple{Int, Float64}}}
    sealing::Vector{Int}
    buffer::Vector{Tuple{Int, Float64}}
    satnum::Vector{Int}
end

mutable struct CSP11ReportHandler
    time::Float64
    next_time::Float64
end

function reformat_weights(weights::AbstractArray{<:Real})
    return [(i, Float64(w)) for (i, w) in enumerate(weights) if w > 0.0]
end

function CSP11ReportHelper(domain; specase = get(domain, :spe11_case, :b))
    points = domain[:spe11_observation_cells]
    buffer_fraction = domain[:spe11_buffer_volumes]./domain[:volumes]
    return CSP11ReportHelper(specase, points,
        (reformat_weights(domain[:A]), reformat_weights(domain[:B])),
        findall(==(1), domain[:satnum]), reformat_weights(buffer_fraction), domain[:satnum])
end

reservoir_reporting_state(state) = haskey(state, :Reservoir) ? state[:Reservoir] : state

"""
    get_reporting_hook(path, domain; time_offset, interval)

Write sparse measurables at the statement's reporting times, including time zero.
The hook reads the current converged state and linearly interpolates measurables
between accepted ministeps when a reporting time lies between them. Failed steps
produce no rows. Convection is intentionally left as NaN.
"""
function get_reporting_hook(path, domain; specase = get(domain, :spe11_case, :b),
        time_offset = get(domain, :spe11_initialization_time, spe11_case_specification(specase).initialization_time),
        interval = spe11_case_specification(specase).sparse_interval)
    isfinite(interval) && interval > 0.0 || throw(ArgumentError("Reporting interval must be finite and positive"))
    mkpath(path)
    file = joinpath(path, "spe11$(specase)_time_series.csv")
    helper = CSP11ReportHelper(domain; specase = specase)
    header = "# t [s], p1 [Pa], p2 [Pa], mobA [kg], immA [kg], dissA [kg], sealA [kg], mobB [kg], immB [kg], dissB [kg], sealB [kg], M_C [m^2], sealTot [kg]"
    open(file, "w") do io
        println(io, specase == :a ? header : header * ", boundTot [kg]")
    end
    handler = CSP11ReportHandler(0.0, 0.0)
    final_time = spe11_case_specification(specase).simulation_time
    function hook(done, report, sim, dt, forces, max_iter, cfg)
        report[:success] || return (done, report)
        previous_time = handler.time - time_offset
        handler.time += report[:dt]
        current_time = handler.time - time_offset
        tolerance = 64*eps(max(handler.time, abs(time_offset), final_time, 1.0))
        if current_time + tolerance >= handler.next_time
            previous = reporting_values(helper, reservoir_reporting_state(sim.storage.state0))
            current = reporting_values(helper, reservoir_reporting_state(sim.storage.state))
            while handler.next_time <= min(current_time, final_time) + tolerance
                fraction = clamp((handler.next_time - previous_time)/(current_time - previous_time), 0.0, 1.0)
                values = previous .+ fraction.*(current .- previous)
                open(file, "a") do io
                    write_reporting_values(io, (handler.next_time, values...))
                end
                handler.next_time += interval
            end
        end
        return (done, report)
    end
    return hook, file
end

function reporting_values(helper::CSP11ReportHelper, state)
    p1, p2 = helper.points
    a = co2_measures(state, helper.regions[1], helper.satnum)
    b = co2_measures(state, helper.regions[2], helper.satnum)
    values = [Jutul.value(state[:Pressure][p1]), Jutul.value(state[:Pressure][p2]),
        a..., b..., NaN, sealing_co2(state, helper.sealing)]
    if helper.specase != :a
        push!(values, sum(w*Jutul.value(state[:TotalMasses][2, c]) for (c, w) in helper.buffer; init = 0.0))
    end
    return values
end

function write_reporting_values(io, values)
    println(io, join((x isa Real ? Jutul.@sprintf("%1.15e", Jutul.value(x)) : string(x) for x in values), ", "))
end

function co2_measures(state, cells, satnum)
    kr, rho = state[:RelativePermeabilities], state[:PhaseMassDensities]
    x, y = state[:LiquidMassFractions], state[:VaporMassFractions]
    saturation, volume = state[:Saturations], state[:FluidVolume]
    mobile = immobile = dissolved = seal = 0.0
    for (cell, weight) in cells
        free = Jutul.value(volume[cell]*saturation[2, cell]*rho[2, cell]*y[2, cell])
        dissolved += weight*Jutul.value(volume[cell]*saturation[1, cell]*rho[1, cell]*x[2, cell])
        if Jutul.value(kr[2, cell]) > 0.0
            mobile += weight*free
        else
            immobile += weight*free
        end
        if satnum[cell] == 1
            seal += weight*Jutul.value(state[:TotalMasses][2, cell])
        end
    end
    return mobile, immobile, dissolved, seal
end

sealing_co2(state, cells) = sum(Jutul.value(state[:TotalMasses][2, c]) for c in cells; init = 0.0)

"""Validate mapping weights and cache normalization for intensive/extensive fields."""
function reporting_mapping(weights, nc, dims)
    all(n -> n isa Integer && n > 0, dims) || throw(ArgumentError("Reporting dimensions must be positive integers"))
    entries = Tuple{Int, Int, Float64}[]
    source_totals, target_totals = zeros(nc), zeros(prod(dims))
    for (i, j, weight) in weights
        i isa Integer && 1 <= i <= length(target_totals) || throw(ArgumentError("Invalid reporting-cell index $i"))
        j isa Integer && 1 <= j <= nc || throw(ArgumentError("Invalid simulation-cell index $j"))
        isfinite(weight) && weight >= 0.0 || throw(ArgumentError("Mapping weights must be finite and non-negative"))
        iszero(weight) && continue
        push!(entries, (i, j, weight))
        source_totals[j] += weight
        target_totals[i] += weight
    end
    all(>(0.0), source_totals) || throw(ArgumentError("Every simulation cell must be covered by the reporting map"))
    return (entries = entries, source_totals = source_totals, target_totals = target_totals, dims = dims)
end

function map_to_reporting_grid(state, mapping; temperature = nothing)
    state = reservoir_reporting_state(state)
    if haskey(state, :Temperature)
        temperature = state[:Temperature]
    elseif isnothing(temperature)
        throw(ArgumentError("Temperature must be supplied for an isothermal state that omits it"))
    end
    fields = (state[:Pressure], view(state[:Saturations], 2, :),
        view(state[:LiquidMassFractions], 2, :), view(state[:VaporMassFractions], 1, :),
        view(state[:PhaseMassDensities], 1, :), view(state[:PhaseMassDensities], 2, :),
        temperature, view(state[:TotalMasses], 2, :))
    mapped = ntuple(_ -> zeros(length(mapping.target_totals)), 8)
    for (i, j, weight) in mapping.entries
        for field in 1:7
            mapped[field][i] += weight/mapping.target_totals[i]*Jutul.value(fields[field][j])
        end
        mapped[8][i] += weight/mapping.source_totals[j]*Jutul.value(fields[8][j])
    end
    for i in eachindex(mapping.target_totals)
        if iszero(mapping.target_totals[i])
            for field in 1:7
                mapped[field][i] = NaN
            end
        end
    end
    return Dict(zip((:p, :sg, :x_co2, :y_h2o, :rho_w, :rho_g, :T, :co2_mass), mapped))
end

function map_to_reporting_grid(case, states::AbstractVector; weights, dims)
    mapping = reporting_mapping(weights, number_of_cells(case.model[:Reservoir].domain), dims)
    temperature = hasproperty(case, :parameters) ? get(case.parameters[:Reservoir], :Temperature, nothing) : nothing
    return [map_to_reporting_grid(state, mapping; temperature = temperature) for state in states]
end

"""
    write_reporting_grid(case, states, path, specase; weights, times, dims)

Map all cells conservatively and write the statement's dense reporting times.
`dims` uses (x, y, elevation), with x varying fastest in the mapping indices.
`times`, when supplied, is measured in seconds since injection. Include the
injection-time state; interpolation between supplied snapshots is linear.
"""
function write_reporting_grid(case, states, path, specase::Symbol; weights,
        dims = spe11_case_specification(specase).report_dims, times = nothing,
        time_offset = get(reservoir_domain(case.model), :spe11_initialization_time,
            spe11_case_specification(specase).initialization_time))
    spec = spe11_case_specification(specase)
    if isnothing(times)
        times = cumsum(case.dt) .- time_offset
        if length(states) == length(times) + 1
            times = [-time_offset; times]
        end
    end
    length(states) == length(times) || throw(ArgumentError("There must be one time per state"))
    isempty(times) && throw(ArgumentError("At least one state is required"))
    all(isfinite, times) && all(>(0.0), diff(times)) || throw(ArgumentError("Snapshot times must be finite and strictly increasing"))
    tolerance = 64*eps(max(maximum(abs, times), spec.simulation_time, 1.0))
    first(times) <= tolerance && last(times) >= -tolerance || throw(ArgumentError("Dense reporting requires the state at injection time zero"))
    mapped = map_to_reporting_grid(case, states; weights = weights, dims = dims)
    mkpath(path)
    for time in spe11_reporting_times(specase)
        time > last(times) + tolerance && break
        index = searchsortedlast(times, time)
        if index == 0 && abs(first(times) - time) <= tolerance
            data = first(mapped)
        elseif index == length(times) || abs(times[index] - time) <= tolerance
            data = mapped[index]
        else
            fraction = (time - times[index])/(times[index+1] - times[index])
            data = Dict(key => mapped[index][key] .+ fraction.*(mapped[index+1][key] .- mapped[index][key]) for key in keys(mapped[index]))
        end
        unit = specase == :a ? si_unit(:hour) : spe11_year
        stamp, suffix = round(Int, time/unit), specase == :a ? "h" : "y"
        file = joinpath(path, "spe11$(specase)_spatial_map_$(stamp)$(suffix).csv")
        write_timestep_dense(file, dims, spec.physical_size, data;
            thermal = specase != :a, three_dimensional = specase == :c)
    end
end

function write_timestep_dense(file, dims, physical_size, data; thermal, three_dimensional)
    nx, ny, nz = dims
    dx, dy, dz = physical_size./dims
    open(file, "w") do io
        header = three_dimensional ? "# x [m], y [m], z [m]" : "# x [m], z [m]"
        header *= ", pressure [Pa], gas saturation [-], mass fraction of CO2 in liquid [-], mass fraction of H2O in vapor [-], phase mass density gas [kg/m3], phase mass density water [kg/m3], total mass CO2 [kg]"
        println(io, thermal ? header * ", temperature [C]" : header)
        index = 1
        for k in 1:nz, j in 1:ny, i in 1:nx
            coordinates = three_dimensional ? ((i - 0.5)*dx, (j - 0.5)*dy, (k - 0.5)*dz) : ((i - 0.5)*dx, (k - 0.5)*dz)
            values = (coordinates..., (data[key][index] for key in (:p, :sg, :x_co2, :y_h2o, :rho_g, :rho_w, :co2_mass))...)
            write_reporting_values(io, thermal ? (values..., data[:T][index] - 273.15) : values)
            index += 1
        end
    end
end

function print_performance(path, states, reports, specase::Symbol; model = nothing,
        time_offset = spe11_case_specification(specase).initialization_time)
    mkpath(path)
    degrees_of_freedom = isnothing(model) ? "n/a" : Jutul.number_of_degrees_of_freedom(model)
    current_time = total_time = 0.0
    open(joinpath(path, "spe11$(specase)_performance_time_series_detailed.csv"), "w") do io
        println(io, "# t [s], avg timestep [s], failed timesteps, mass of CO2 [kg], ndof, nliter, nres, liniter, runtime [s], tlinsol [s]")
        for (i, report) in enumerate(reports)
            current_time += report_timesteps([report])[1]
            total_time += report[:total_time]
            current_time < time_offset && continue
            ministeps = report[:ministeps]
            successful = filter(step -> step[:success], ministeps)
            average_step = isempty(successful) ? "n/a" : sum(step[:dt] for step in successful)/length(successful)
            stats = [step[:stats] for step in ministeps]
            state = reservoir_reporting_state(states[i])
            values = (current_time - time_offset, average_step, count(step -> !step[:success], ministeps),
                sum(state[:TotalMasses][2, :]), degrees_of_freedom,
                sum(s.newtons for s in stats), sum(s.linearizations for s in stats),
                sum(s.linear_iterations for s in stats), total_time,
                sum(s.linear_setup + s.linear_solve for s in stats))
            write_reporting_values(io, values)
        end
    end
end
