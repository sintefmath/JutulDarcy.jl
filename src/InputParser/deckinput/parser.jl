include("units.jl")
include("utils.jl")
include("keywords/keywords.jl")

function parse_deck_file(filename; kwarg...)
    outer_data = Dict{String, Any}()
    parse_deck_file!(outer_data, filename; kwarg...)
    return outer_data
end

function parse_deck_file!(outer_data, filename, data = outer_data;
        skip_mode = false,
        verbose = false,
        sections = [:RUNSPEC, :GRID, :PROPS, :REGIONS, :SOLUTION, :SCHEDULE],
        skip = [:SUMMARY],
        units::Union{Symbol, Nothing} = :si,
        input_units::Union{Symbol, Nothing} = nothing,
        unit_systems = missing
    )
    function get_unit_system_pair(from::Symbol, target::Symbol)
        from_sys = DeckUnitSystem(from)
        target_sys = DeckUnitSystem(target)
        # We have a pair of unit systems to convert between
        return (from = from_sys, to = target_sys)
    end

    if isnothing(units)
        # nothing for units means no conversion
        @assert ismissing(unit_systems)
        unit_systems = nothing
    end
    if !isnothing(units) && !isnothing(input_units)
        # Input and output units are both specified, potentially overriding
        # what's in RUNSPEC. This will also check that they are valid symbols
        # for us.
        unit_systems = get_unit_system_pair(input_units, units)
    end
    filename = realpath(filename)
    basedir = dirname(filename)
    function msg(s)
        if verbose
            println("PARSER: $s")
        end
    end
    msg("Starting parse of $filename...")
    f = open(filename, "r")

    try
        allsections = vcat(sections, skip)
        while !eof(f)
            m = next_keyword!(f)
            if isnothing(m)
                # Do nothing
            elseif m in allsections
                msg("Starting section $m")
                data = new_section(outer_data, m)
                skip_mode = m in skip
                # Check if we have passed RUNSPEC so that units can be set
                runspec_passed = m != :RUNSPEC && haskey(outer_data, "RUNSPEC")
                unit_system_not_initialized = ismissing(unit_systems)
                if runspec_passed && unit_system_not_initialized
                    unit_systems = get_unit_system_pair(current_unit_system(outer_data), units)
                end
            elseif m == :INCLUDE
                next = strip(readline(f))
                if endswith(next, '/')
                    next = rstrip(next, '/')
                else
                    readline(f)
                end
                include_path = clean_include_path(basedir, next)
                msg("Including file $include_path...")
                parse_deck_file!(
                    outer_data, include_path, data,
                    verbose = verbose,
                    sections = sections,
                    skip = skip,
                    skip_mode = skip_mode,
                    units = units,
                    unit_systems = unit_systems
                )
            elseif m == :END
                # All done!
                break
            elseif !skip_mode
                parse_keyword!(data, outer_data, unit_systems, f, Val(m))
            end
        end
    finally
        close(f)
    end
    return outer_data
end

