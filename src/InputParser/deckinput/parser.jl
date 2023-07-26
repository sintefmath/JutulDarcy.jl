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
        unit_systems = missing
    )
    if isnothing(units)
        # nothing for units means no conversion
        @assert ismissing(unit_systems)
        unit_systems = nothing
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
                from_sys = DeckUnitSystem(current_unit_system(outer_data))
                target_sys = DeckUnitSystem(units)
                # We have a pair of unit systems to convert between
                unit_systems = (from = from_sys, to = target_sys)
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
    close(f)
    # if !isnothing(units)
    #     convert_deck_units!(outer_data, units)
    # end
    return outer_data
end

