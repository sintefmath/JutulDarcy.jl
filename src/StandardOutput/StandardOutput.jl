module StandardOutput

using Dates: DateTime
using Jutul: JutulCase, UnstructuredMesh, evaluate_all_secondary_variables,
    number_of_cells, physical_representation, report_times
using GeoEnergyIO: BinaryKeyword, KeywordFile, write_egrid, write_restart, write_summary
using GeoEnergyIO.InputParser: DeckUnitSystem
using GeoEnergyIO.StandardIO: case_base, convert_array, unit_code
using GeoEnergyIO.StandardIO.output: fake_grid, grid_dims, restart_block, standard_headers
using MultiComponentFlash: molar_masses
using ..JutulDarcy: ReservoirSimResult, AqueousPhase, LiquidPhase, VaporPhase,
    CompositionalSystem, deck_output_units, get_phases, get_reference_phase_index,
    mass_to_mole_fractions, reservoir_domain, reservoir_model, summary_result

export write_standard_output, write_standard_summary, write_standard_restart, write_standard_egrid

# Eclipse PRESSURE uses oil, gas or water pressure in this order of preference.
const RESTART_PHASES = (
    (phase = LiquidPhase(), saturation = "SOIL", mask = 1, density = "OIL_DEN", viscosity = "OIL_VISC", kr = "OILKR", invb = "1OVERBO"),
    (phase = VaporPhase(), saturation = "SGAS", mask = 4, density = "GAS_DEN", viscosity = "GAS_VISC", kr = "GASKR", invb = "1OVERBG"),
    (phase = AqueousPhase(), saturation = "SWAT", mask = 2, density = "WAT_DEN", viscosity = "WAT_VISC", kr = "WATKR", invb = "1OVERBW"),
)

"""
    write_standard_output(result, case, path; summary=true, restart=true, egrid=true,
                          units=missing, unified=true)

Write industry-standard summary, reservoir restart and EGRID files for a
[`ReservoirSimResult`](@ref) and its `JutulCase`. `path` is the common output
basename, including an optional directory. Set any of `summary`, `restart` or
`egrid` to `false` to skip it. Set `unified=false` to write separate summary
and restart data files.

The deck's `RUNSPEC` selects the output units when available, with metric as
the default. Use `units=:metric`, `:field` or `:lab` to override that choice.
The result is a named tuple of the written paths, with `nothing` for any
skipped output.
"""
function write_standard_output(result::ReservoirSimResult, case::JutulCase, path::AbstractString;
        summary::Bool = true, restart::Bool = true, egrid::Bool = true,
        units = missing, unified::Bool = true)
    summary_paths = summary ? write_standard_summary(result, case, path; units, unified) : nothing
    restart_paths = restart ? write_standard_restart(result, case, path; units, unified) : nothing
    egrid_path = egrid ? write_standard_egrid(case, path; units) : nothing
    return (summary = summary_paths, restart = restart_paths, egrid = egrid_path)
end

"""
    write_standard_summary(result, case, path; units=missing, unified=true)

Write SMSPEC and UNSMRY (or separate Snnnn files) from the field, well and
group summary vectors. Return the SMSPEC path and the summary data paths.
"""
function write_standard_summary(result::ReservoirSimResult, case::JutulCase, path::AbstractString;
        units = missing, unified::Bool = true)
    units = output_units(case, units)
    layout = grid_layout(case)
    model = reservoir_model(case)
    states = result.states
    if haskey(model.secondary_variables, :Saturations) && any(s -> !haskey(s, :Saturations), states)
        parameters = case.parameters[:Reservoir]
        states = map(s -> complete_output_state(model, s, parameters, (:Saturations,)), states)
    end
    summary = summary_result(case.model, result.wells, states, units;
        input_data = case.input_data, start_date = output_start_date(case))
    summary["DIMENS"] = layout.dims
    return write_summary(summary, path; unified)
end

"""
    write_standard_egrid(case, path; units=missing)

Write and return the EGRID path. Use the geometry in `case.input_data["GRID"]`
when available. Otherwise, write a synthetic `nc × 1 × 1` grid with one unit
cube per reservoir cell, in simulation order.
"""
function write_standard_egrid(case::JutulCase, path::AbstractString; units = missing)
    units = output_units(case, units)
    (; grid, dims, cell_map) = grid_layout(case)
    if isnothing(grid)
        grid = fake_grid(dims[1])
    else
        grid = Dict{String, Any}(grid)
        mask = falses(dims)
        mask[cell_map] .= true
        grid["ACTNUM"] = mask
    end
    return write_egrid(grid, case_base(path); units)
end

"""
    write_standard_restart(result, case, path; units=missing, unified=true)

Write UNRST (or separate Xnnnn files) with reservoir solution arrays for each
completed report step. Include pressure, phase saturations, active black-oil
`RS` and `RV`, and temperature when defined. Compositional output includes
numbered `ZMF`, `XMF` and `YMF` mole fractions in equation-of-state component
order. Immiscible water is represented separately by `SWAT`. Phase densities,
viscosities, relative permeabilities and inverse formation volume factors are
included when those properties were saved in the simulation result.
Select them with `extra_outputs` when setting up the reservoir model.

The files contain reservoir solution data for post-processing. Return the
unified file path or the separate paths.
"""
function write_standard_restart(result::ReservoirSimResult, case::JutulCase, path::AbstractString;
        units = missing, unified::Bool = true)
    units = output_units(case, units)
    (; dims, cell_map) = grid_layout(case)
    # Restart arrays are stored in the order of the active cells in ACTNUM.
    order = issorted(cell_map) ? eachindex(cell_map) : sortperm(cell_map)
    model = reservoir_model(case)
    system = model.system
    phases = get_phases(system)
    parameters = case.parameters[:Reservoir]
    deck_units = DeckUnitSystem(units)
    active_phases = [(findfirst(isequal(spec.phase), phases), spec) for spec in RESTART_PHASES if spec.phase in phases]
    phase_mask = sum(p -> last(p).mask, active_phases; init = 0)
    ref_phase = get_reference_phase_index(system)
    output_phase = isempty(active_phases) ? ref_phase : first(active_phases)[1]
    core = (:Saturations, :Rs, :Rv, :LiquidMassFractions, :VaporMassFractions)
    if output_phase != ref_phase
        core = (core..., :CapillaryPressure)
    end
    required = filter(key -> haskey(model.secondary_variables, key), core)
    masses = system isa CompositionalSystem ? molar_masses(system.equation_of_state) : nothing
    start_date = output_start_date(case)

    states = result.result.states
    isempty(states) && throw(ArgumentError("Restart output requires at least one completed report step."))
    elapsed = report_times(result.result.reports[eachindex(states)])
    steps = RestartSteps(length(states)) do i
        stored = states[i][:Reservoir]
        state = complete_output_state(model, stored, parameters, required)
        headers = standard_headers(dims, length(cell_map);
            units, phases = phase_mask, elapsed = elapsed[i], start_date, report_step = i)
        solution = BinaryKeyword[]
        pressure = view(state[:Pressure], order)
        if output_phase != ref_phase && haskey(state, :CapillaryPressure)
            row = output_phase - (output_phase > ref_phase)
            pressure = pressure .+ view(state[:CapillaryPressure], row, order)
        end
        push!(solution, BinaryKeyword("PRESSURE", Float32.(convert_array("PRESSURE", pressure, units; to_si = false))))
        for (key, name) in ((:Rs, "RS"), (:Rv, "RV"))
            if haskey(state, key)
                values = convert_array(name, view(state[key], order), units; to_si = false)
                push!(solution, BinaryKeyword(name, Float32.(values)))
            end
        end
        temperature = get(state, :Temperature, get(parameters, :Temperature, nothing))
        if !isnothing(temperature)
            values = convert_array("TEMP", view(temperature, order), units; to_si = false)
            push!(solution, BinaryKeyword("TEMP", Float32.(values)))
        end
        for (index, spec) in active_phases
            (; saturation, density, viscosity, kr, invb) = spec
            if haskey(state, :Saturations)
                push!(solution, BinaryKeyword(saturation, Float32.(view(state[:Saturations], index, order))))
            elseif length(phases) == 1
                push!(solution, BinaryKeyword(saturation, ones(Float32, length(cell_map))))
            end
            fvf = spec.phase isa VaporPhase ? deck_units.gas_formation_volume_factor : deck_units.liquid_formation_volume_factor
            for (key, keyword, scale) in ((:PhaseMassDensities, density, inv(deck_units.density)),
                    (:PhaseViscosities, viscosity, inv(deck_units.viscosity)),
                    (:RelativePermeabilities, kr, 1.0), (:ShrinkageFactors, invb, fvf))
                if haskey(stored, key)
                    push!(solution, BinaryKeyword(keyword, Float32.(view(stored[key], index, order).*scale)))
                end
            end
        end
        if !isnothing(masses)
            append_mole_fractions!(solution, state, order, masses)
        end
        return restart_block(headers, solution; report_step = i)
    end
    return write_restart(steps, path; unified)
end

function complete_output_state(model, state, parameters, required)
    missing_keys = filter(key -> !haskey(state, key), required)
    isempty(missing_keys) && return state
    report_parameters = Dict(key => get(state, key, v) for (key, v) in pairs(parameters))
    evaluated = evaluate_all_secondary_variables(model, state, report_parameters)
    completed = copy(state)
    for key in missing_keys
        completed[key] = evaluated[key]
    end
    return completed
end

function append_mole_fractions!(solution, state, order, masses)
    overall = state[:OverallMoleFractions]
    for i in eachindex(masses)
        push!(solution, BinaryKeyword("ZMF$i", Float32.(view(overall, i, order))))
    end
    ncomponents = Val(length(masses))
    for (key, prefix) in ((:LiquidMassFractions, "XMF"), (:VaporMassFractions, "YMF"))
        fractions = state[key]
        mole = [mass_to_mole_fractions(fractions, masses, c, ncomponents) for c in order]
        for i in eachindex(masses)
            push!(solution, BinaryKeyword("$prefix$i", Float32[x[i] for x in mole]))
        end
    end
end

# GeoEnergyIO takes a vector of report blocks. Build each block as it is written.
struct RestartSteps{F} <: AbstractVector{KeywordFile}
    make_step::F
    count::Int
end

Base.size(steps::RestartSteps) = (steps.count,)
Base.IndexStyle(::Type{<:RestartSteps}) = IndexLinear()
function Base.getindex(steps::RestartSteps, i::Int)
    @boundscheck checkbounds(steps, i)
    return steps.make_step(i)
end

function input_section(case, name)
    data = case.input_data
    return data isa AbstractDict ? get(data, name, nothing) : nothing
end

function output_units(case, units)
    if ismissing(units)
        units = deck_output_units(input_section(case, "RUNSPEC"))
    end
    units = Symbol(lowercase(string(units)))
    unit_code(units)
    return units
end

function output_start_date(case)
    if !isnothing(case.start_date)
        return case.start_date
    end
    runspec = input_section(case, "RUNSPEC")
    return isnothing(runspec) ? DateTime(1970) : get(runspec, "START", DateTime(1970))
end

function grid_layout(case)
    mesh = physical_representation(reservoir_domain(case))
    nc = number_of_cells(mesh)
    grid = input_section(case, "GRID")
    if isnothing(grid)
        return (; grid, dims = (nc, 1, 1), cell_map = Base.OneTo(nc))
    end
    dims = grid_dims(grid)
    if mesh isa UnstructuredMesh && mesh.structure isa CartesianIndex
        cell_map = mesh.cell_map
    else
        actnum = get(grid, "ACTNUM", nothing)
        cell_map = isnothing(actnum) ? Base.OneTo(prod(dims)) : findall(!iszero, vec(actnum))
    end
    length(cell_map) == nc || throw(ArgumentError("The input GRID does not match the reservoir's active cells."))
    return (; grid, dims, cell_map)
end

end
