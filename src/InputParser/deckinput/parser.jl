include("units.jl")
include("utils.jl")
include("keywords/keywords.jl")

Base.@kwdef struct ParserVerbosityConfig
    silent::Bool = false
    verbose::Bool = true
    warn_feature::Bool = true
    warn_parsing::Bool = true
    warn_limit::Int = 1
    warn_count::Dict{String, Int} = Dict{String, Int}()
end

@enum PARSER_WARNING PARSER_MISSING_SUPPORT PARSER_JUTULDARCY_MISSING_SUPPORT PARSER_JUTULDARCY_PARTIAL_SUPPORT PARSER_PARTIAL_SUPPORT

function keyword_header(current_data, keyword)
    if haskey(current_data, "CURRENT_SECTION")
        section = current_data["CURRENT_SECTION"]
        out = "$section/$keyword"
    else
        out = "$keyword"
    end
    return out
end

function parser_message(cfg::ParserVerbosityConfig, outer_data, keyword, msg::String; important = false)
    if !cfg.silent
        if important || cfg.verbose
            println("$(keyword_header(outer_data, keyword)): $msg")
        end
    end
end

function parser_message(cfg::ParserVerbosityConfig, outer_data, keyword, msg::PARSER_WARNING, reason = missing)
    if cfg.silent
        return
    end
    if msg == PARSER_MISSING_SUPPORT
        text_msg = "Unsupported keyword. It will be ignored."
        do_print = cfg.warn_parsing
    elseif msg == PARSER_JUTULDARCY_MISSING_SUPPORT
        text_msg = "$keyword is not supported by JutulDarcy solvers. It will be ignored in simulations."
        do_print = cfg.warn_feature
    elseif msg == PARSER_JUTULDARCY_PARTIAL_SUPPORT
        text_msg = "$keyword is only partially supported by JutulDarcy solvers."
        do_print = cfg.warn_feature
    else
        @assert msg == PARSER_PARTIAL_SUPPORT
        text_msg = "Parser has only partial support. $keyword may have missing or wrong entries."
        do_print = cfg.warn_parsing
    end
    if !ismissing(reason)
        text_msg = "$text_msg\n\n$reason"
    end
    if haskey(cfg.warn_count, keyword)
        cfg.warn_count[keyword] += 1
    else
        cfg.warn_count[keyword] = 1
    end
    if cfg.warn_count[keyword] <= cfg.warn_limit
        jutul_message("Parser", "$(keyword_header(outer_data, keyword)) - $text_msg", color = :yellow)
    end
end

"""
    parse_data_file(filename; unit = :si)
    data = parse_data_file("MY_MODEL.DATA")

Parse a .DATA file given by the `filename` (industry standard input file) into a
Dict. Units will be converted to strict SI unless you pass something else like
`units = :field`. Setting `units = nothing` will skip unit conversion. Note that
JutulDarcy assumes that the unit system is internally consistent. It is highly
recommended to parse to the SI units if you want to perform simulations.

The best publicly available documentation on this format is available from the
Open Porous Media (OPM) project's webpages: [OPM Flow manual
](https://opm-project.org/?page_id=955).

# Keyword arguments
- `warn_parsing=true`: Produce a warning when keywords are not supported (or
  partially supported) by the parser.
- `warn_feature`=true`: Produce a warning when keywords are supported, but have
  limited or missing support in the numerical solvers.
- `units=:si`: Symbol that indicates the unit system to be used in the output.
  Setting this to `nothing` will return values without conversion, i.e. exactly
  what is in the input files. `:si` will use strict SI. Other alternatives are
  `:field` and `:metric`. `:lab` is currently unsupported.

# Note
This function is experimental and only covers a small portion of the keywords
that exist for various simulators. You will get warnings that indicate the level
of support for keywords in both the parser and the numerical solvers when known
keywords with limited support. Pull requests for new keywords are welcome!
"""
function parse_data_file(filename; kwarg...)
    outer_data = Dict{String, Any}()
    parse_data_file!(outer_data, filename; kwarg...)
    delete!(outer_data, "CURRENT_SECTION")
    return outer_data
end

function parse_data_file!(outer_data, filename, data = outer_data;
        skip_mode = false,
        verbose = false,
        sections = [:RUNSPEC, :GRID, :PROPS, :REGIONS, :SOLUTION, :SCHEDULE, :EDIT],
        skip = [:SUMMARY],
        units::Union{Symbol, Nothing} = :si,
        warn_parsing = true,
        warn_feature = true,
        basedir = missing,
        silent = false,
        is_outer = true,
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
    if ismissing(basedir)
        basedir = dirname(filename)
    end
    f = open(filename, "r")

    cfg = ParserVerbosityConfig(
        verbose = verbose,
        warn_parsing = warn_parsing,
        warn_feature = warn_feature,
        silent = silent
    )
    parser_message(cfg, outer_data, "PARSER", "Starting parse of $filename...")
    try
        allsections = vcat(sections, skip)
        while !eof(f)
            m = next_keyword!(f)
            if isnothing(m)
                continue
            end
            parsed_keyword = false
            parser_message(cfg, outer_data, "$m", "Starting parse...")
            t_p = @elapsed if m == :SKIP || m == :ENDSKIP
                if m == :SKIP
                    skip_mode = true
                else
                    @assert skip_mode
                    skip_mode = false
                end
            elseif m in allsections
                # Check if we have passed RUNSPEC so that units can be set
                runspec_passed = m != :RUNSPEC && haskey(outer_data, "RUNSPEC")
                unit_system_not_initialized = ismissing(unit_systems)
                if runspec_passed && unit_system_not_initialized
                    unit_systems = get_unit_system_pair(current_unit_system(outer_data), units)
                end
                finish_current_section!(data, unit_systems, cfg, outer_data)
                parser_message(cfg, outer_data, "$m", "Starting new section.")
                data = new_section(outer_data, m)
                skip_mode = m in skip
            elseif m == :INCLUDE
                next = strip(readline(f))
                if endswith(next, '/')
                    next = rstrip(next, '/')
                else
                    readline(f)
                end
                include_path = clean_include_path(basedir, next)
                parser_message(cfg, outer_data, "$m", "Including file: $include_path. Basedir: $basedir with INCLUDE = $next)")
                parse_data_file!(
                    outer_data, include_path, data,
                    verbose = verbose,
                    warn_parsing = warn_parsing,
                    warn_feature = warn_feature,
                    silent = silent,
                    sections = sections,
                    skip = skip,
                    basedir = basedir,
                    skip_mode = skip_mode,
                    is_outer = false,
                    units = units,
                    unit_systems = unit_systems
                )
            elseif m in (:DATES, :TIME, :TSTEP)
                parsed_keyword = true
                parse_keyword!(data, outer_data, unit_systems, cfg, f, Val(m))
                # New control step starts after this
                data = OrderedDict{String, Any}()
                push!(outer_data["SCHEDULE"]["STEPS"], data)
            elseif m == :END
                # All done!
                finish_current_section!(data, unit_systems, cfg, outer_data)
                break
            elseif skip_mode
                parser_message(cfg, outer_data, "$m", "Keyword skipped.")
            else
                parse_keyword!(data, outer_data, unit_systems, cfg, f, Val(m))
                parsed_keyword = true
            end
            if parsed_keyword
                parser_message(cfg, outer_data, "$m", "Keyword parsed in $(t_p)s.")
            end
        end
    finally
        if is_outer
            finish_current_section!(data, unit_systems, cfg, outer_data)
        end
        close(f)
    end
    return outer_data
end

"""
    parse_grdecl_file("mygrid.grdecl"; actnum_path = missing, kwarg...)

Parse a GRDECL file separately from the full input file. Note that the GRID
section does not contain units - passing the `input_units` keyword is therefore
highly recommended.

# Keyword arguments
 - `units=:si`: Units to use for return values. Requires `input_units` to be set.
 - `input_units=nothing`: The units the file is given in.
 - `verbose=false`: Toggle verbosity.

"""
function parse_grdecl_file(filename; actnum_path = missing, kwarg...)
    outer_data = Dict{String, Any}()
    data = new_section(outer_data, :GRID)
    parse_data_file!(outer_data, filename, data; kwarg...)
    if !ismissing(actnum_path)
        parse_data_file!(outer_data, actnum_path, data; kwarg...)
    end
    if !haskey(data, "ACTNUM")
        data["ACTNUM"] = fill(true, data["cartDims"])
    end
    delete!(data, "CURRENT_BOX")
    return data
end

function finish_current_section!(data, units, cfg, outer_data)
    if haskey(outer_data, "CURRENT_SECTION")
        v = outer_data["CURRENT_SECTION"]
        v::Symbol
        finish_current_section!(data, units, cfg, outer_data, Val(v))
    end
end
