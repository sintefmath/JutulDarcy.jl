using Base64: base64decode

function dir_to_doc_name(x::String)
    x = replace(x, "_" => " ")
    x = uppercase(x[1:1])*x[2:end]
    return x
end

function get_example_paths(; check_empty = true)
    basepth = joinpath(@__DIR__, "..", "examples")
    examples = OrderedDict()
    examples["introduction"] = []
    examples["workflow"] = []
    examples["data_assimilation"] = []
    examples["geothermal"] = []
    examples["compositional"] = []
    examples["fractures"] = []
    examples["discretization"] = []
    examples["properties"] = []
    examples["validation"] = []
    for excat in readdir(basepth)
        if isdir(joinpath(basepth, excat))
            for exfile in readdir(joinpath(basepth, excat))
                if endswith(exfile, ".jl")
                    if !haskey(examples, excat)
                        examples[excat] = []
                    end
                    filename = first(splitext(exfile))
                    push!(examples[excat], filename)
                end
            end
        end
    end
    if check_empty
        for (k, v) in pairs(examples)
            @assert length(examples[k]) > 0 "No examples found for category $k"
        end
    end
    return examples
end

function example_path_jl(cname, pth)
    # Base directory
    jutul_dir = realpath(joinpath(@__DIR__, ".."))
    return joinpath(jutul_dir, "examples", cname, "$pth.jl")
end

function parse_tags(text)
    # Match the pattern <tags: ...> and extract the content
    m = match(r"<tags:\s*([^>]+)>", text)
    if !isnothing(m)
        # Split by comma and strip whitespace
        tags = strip.(split(m.captures[1], ","))
        return tags
    end
    return nothing
end

function example_tags(pth)
    tags = String[]
    lines = readlines(pth)
    for line in lines
        t = parse_tags(line)
        if !isnothing(t)
            append!(tags, t)
            break
        end
    end
    return tags
end


# exlist = get_example_paths(check_empty = false)

function all_tags()
    descr = OrderedDict()
    descr["Introduction"] = "Examples that illustrate basic features of JutulDarcy.jl and how to get started with the simulator. These examples are a good place to start because they are more detailed and easier to follow than the other examples where it is assumed that you are familiar with the features used in most examples."
    descr["StartToFinish"] = "Tag for examples that go through model setup in detail from start to finish, including meshing, setting up properties, wells, and running a simulation. These examples are a good place to start if you want to see how to set up a complete model."
    descr["Advanced"] = "Examples that illustrate more advanced features of JutulDarcy.jl. These examples assume that you are already familiar with the basics of the simulator (e.g. how a reservoir is represented, how wells are set up) as little attention is given to the basics."
    descr["Validation"] = "These examples validate the simulator on well-known problems by comparing results to other simulators or analytical solutions."

    descr["Wells"] = "Examples that have a particular focus on wells and are a good place to look if you are building your own model with wells, both in terms of realizing wells in a model and operating them with different types of constraints."
    descr["InputFile"] = "Examples that illustrate how to set up and run simulations based on input files (e.g. Eclipse/.DATA format)."
    descr["Differentiability"] = "Examples that illustrate how to use the differentiable features of JutulDarcy.jl, including sensitivity calculations and gradient-based optimization."
    descr["HistoryMatching"] = "Demonstrations of how to use JutulDarcy.jl for history matching/data assimilation, including gradient-based optimization of model parameters to match observed data."
    descr["Discretizations"] = "Examples that illustrate how to use different discretizations for flow and transport in JutulDarcy.jl, including advanced discretizations such as multipoint flux approximations (MPFA) and high-resolution schemes (WENO)."
    descr["Meshing"] = "These examples cover meshing (e.g. by calling Gmsh or other packages for mesh generation)."
    descr["MachineLearning"] = "Examples that involve machine learning techniques, such as using neural networks to model relative permeability or other properties, or integrating machine learning type models into the simulation workflow."
    descr["Tracers"] = "Examples that illustrate how to use tracers in JutulDarcy.jl, including setting up passive and active tracers, and using tracers for enhanced oil recovery (EOR) simulations."
    descr["ModelReduction"] = "Examples that illustrate how to use model reduction techniques in JutulDarcy.jl, including coarsening/upscaling, and other methods to reduce the computational cost of simulations."
    descr["Optimization"] = "Examples that illustrate how to use JutulDarcy.jl for optimization problems, including gradient-based optimization of well controls and other parameters to maximize an objective function."
    descr["Properties"] = "Examples that illustrate how to set up and use different dynamic properties in JutulDarcy.jl, including PVT properties, relative permeability, hysteresis and capillary pressure."

    descr["Immiscible"] = "Examples that make use of the immiscible/dead-oil model for PVT descriptions."
    descr["CO2"] = "Examples that involve simulation of geological sequestration  of CO2 (carbon storage and sequestration / CCS), or other types of simulations involving CO2."
    descr["Blackoil"] = "Examples that make use of the blackoil model for PVT descriptions."
    descr["Compositional"] = "Examples that make use of the compositional model for PVT descriptions."
    descr["Geothermal"] = "Examples that simulate recovery and/or storage of heat in the subsurface. See also the dedicated [Fimbul.jl](https://sintefmath.github.io/Fimbul.jl/dev/) module for geothermal simulation with JutulDarcy.jl where additional advanced examples are provided."

    out = OrderedDict()
    colors = to_colormap(:Set1_8)
    i = 1
    rgb_html(x, s) = Int(ceil(getfield(x, s)*255))

    for (k, v) in pairs(descr)
        color = colors[mod1(i, length(colors))]
        r = rgb_html(color, :r)
        g = rgb_html(color, :g)
        b = rgb_html(color, :b)

        out[k] = (desc = v, color = "rgb($(r), $(g), $(b))")
        i += 1
    end
    return out
end

function tag_str(tag::AbstractString)
    return tag_str([tag])
end

function tag_str(tagname::AbstractVector)
    tags = all_tags()
    s = "``` @raw html\n"
    for tag in tagname
        info = tags[tag]
        s *= "<ExampleTag text=\"$tag\" color=\"$(info.color)\" />\n"
    end
    s *= "```\n"
    return s
end

function example_tags()
    expths = get_example_paths(check_empty = false)
    out = Dict{String, Any}()
    for key in keys(all_tags())
        out[key] = Tuple{String, String}[]
    end
    if get(ENV, "JUTULDARCY_DOCS_EXAMPLES_SKIP", "0") == "0"
        for (category, example_set) in pairs(expths)
            for exname in example_set
                pth = example_path_jl(category, exname)
                extags = example_tags(pth)
                @assert length(extags) > 0 "Example $exname in $category has no tags"
                for tag in extags
                    @assert haskey(out, tag) "Example $exname in $category has unknown tag $tag"
                    push!(out[tag], (exname, category))
                end
            end
        end
    end
    return out
end

function write_tags()
    tags = all_tags()
    outdir = joinpath(@__DIR__, "src", "examples", "overview")
    mkpath(outdir)
    outpth = joinpath(outdir, "example_overview.md")
    ex_tags = example_tags()
    open(outpth, "w") do io
        println(io, "# Example overview\n")
        println(io, "JutulDarcy.jl comes with a number of examples that illustrate different features of the simulator. The examples are categorized by tags, and you can find the examples with a specific tag below. Examples can have multiple tags.\n")
        for (tag, info) in pairs(tags)
            println(io, "## $tag\n")
            println(io, tag_str(tag))
            println(io, "$(info.desc)\n")
            println(io, "### Examples with the $(lowercase(tag)) tag:\n")
            if length(keys(ex_tags[tag])) == 0
                println(io, "_No examples with this tag yet._\n")
            else
                for (exname, category) in ex_tags[tag]
                    exlink = joinpath("..", "..", "examples", category, "$exname.md")
                    println(io, "1. [$exname]($exlink) (in $category)")
                end
                println(io, "\n")
            end
        end
    end
    println("Wrote tags to $outpth")
end


##
function timer_str()
    start = "example_t_start = time_ns(); # hide\n"
    stop_1 = "\nt_s = (time_ns() - example_t_start) / 1e9 # hide\n"
    stop = stop_1*"println(\"This example took "*raw"$t_s"*" seconds to complete.\") # hide"
    return (start, stop)
end

function post_run_variables_gc()
    # Attempt of post-processing to GC some objects in temporary @example modules
    # https://discourse.julialang.org/t/delete-a-module/62226
    s =  "\nfunction clear_module!(M::Module)        # hide\n"*
    "    for name ∈ names(M, all=true)        # hide\n"*
    "        if !isconst(M, name)             # hide\n"*
    raw"            @eval M $name = $nothing     # hide"*
    "\n        end                              # hide\n"*
    "    end                                  # hide\n"*
    "end                                      # hide\n"*
    "clear_module!(@__MODULE__)               # hide\n"*
    "GC.gc();                                 # hide\n"
    return s
end

function example_info_footer(subdir, exname)
    return "\n\n# ## Example on GitHub\n"*
    "# If you would like to run this example yourself, it can be downloaded from "*
    "the JutulDarcy.jl GitHub repository [as a script](https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/$subdir/$exname.jl)"
    # return "\n\n# ## Example on GitHub\n"*
    # "# If you would like to run this example yourself, it can be downloaded from "*
    # "the JutulDarcy.jl GitHub repository [as a script](https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/$subdir/$exname.jl), "*
    # "or as a [Jupyter Notebook](https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/$subdir/$exname.ipynb)"
end

function update_footer(content, subdir, exname; execute = false)
    info_footer = example_info_footer(subdir, exname)
    gc_footer = post_run_variables_gc()
    start, stop = timer_str()
    # Literate's markdown executor discards its captured explicit displays.
    # Match Documenter's textual display output for the display(x) calls in examples.
    display_setup = execute ? "display(x) = (show(stdout, MIME(\"text/plain\"), x); println()) # hide\n" : ""
    new_content = string(display_setup, start, content, info_footer, stop, gc_footer)
    # Documenter accepts both '# hide' and '#hide'; Literate execution strips only
    # the latter. Normalize comments so hidden setup and cleanup stay hidden.
    if execute
        new_content = replace(new_content, r"#[\t ]+hide[\t \r]*$"m => "#hide")
    end
    return new_content
end

function build_literate_example(category, exname; execute = true)
    in_pth = example_path_jl(category, exname)
    out_dir = joinpath(@__DIR__, "src", "examples", category)
    marker = joinpath(out_dir, "$exname.sha256")
    rm(marker; force = true)
    upd(content) = update_footer(content, category, exname; execute = execute)
    function fixt(content)
        content = replace_tags(content, category, exname)
        if execute
            content = materialize_html_images(content, out_dir, exname)
            content = limit_plain_outputs(content)
        end
        return content
    end
    elapsed = @elapsed Literate.markdown(in_pth, out_dir;
        execute = execute, preprocess = upd, postprocess = fixt,
        edit_commit = "main",
    )
    if execute
        # Only mark successful execution; --docs-only rejects stale pages.
        write(marker, example_fingerprint(category, exname))
    end
    @info "Built example" category exname elapsed
    return nothing
end

function limit_plain_outputs(content; max_lines = 80, max_columns = 160)
    # Literate does not set IOContext(:limit => true) for text/plain results.
    # Track all fences so closing Julia/HTML fences cannot be mistaken for the
    # opening fence of a plain output block.
    output = String[]
    block = String[]
    fence = ""
    limit = false
    for line in split(content, '\n')
        m = match(r"^(`{3,})(.*)$", line)
        if isempty(fence)
            if !isnothing(m)
                fence = m.captures[1]
                limit = fence == "````" && isempty(strip(m.captures[2]))
            end
            push!(output, line)
        elseif !isnothing(m) && length(m.captures[1]) >= length(fence) && isempty(strip(m.captures[2]))
            if limit
                truncated = length(block) > max_lines || any(line -> length(line) > max_columns, block)
                append!(output, truncated ? [first(line, max_columns) for line in Iterators.take(block, max_lines)] : block)
                truncated && push!(output, "… (output truncated)")
                empty!(block)
            end
            push!(output, line)
            fence = ""
            limit = false
        elseif limit
            push!(block, line)
        else
            push!(output, line)
        end
    end
    append!(output, block) # Preserve an unclosed output fence verbatim.
    return join(output, '\n')
end

function materialize_html_images(content, out_dir, exname)
    # Makie's HTML display embeds PNGs as base64. Store these as assets instead
    # of putting megabytes of image data into the Markdown/Vue compilation step.
    index = 0
    pattern = r"src=\"data:image/(png|jpeg);base64,\s*([A-Za-z0-9+/=\s]+)\""
    return replace(content, pattern => function (html)
        m = match(pattern, html)
        index += 1
        name = "$exname-inline-$index.$(m.captures[1])"
        write(joinpath(out_dir, name), base64decode(m.captures[2]))
        return "src=\"./$name\""
    end)
end

function example_fingerprint(category, exname)
    files = [example_path_jl(category, exname), @__FILE__, joinpath(@__DIR__, "make.jl"),
        joinpath(@__DIR__, "example_worker.jl"),
        joinpath(@__DIR__, "Project.toml")]
    manifest = joinpath(@__DIR__, "Manifest.toml")
    isfile(manifest) && push!(files, manifest)
    # Source changes can affect an example without changing the script itself.
    for dir in ("src", "ext")
        for (root, _, names) in walkdir(joinpath(@__DIR__, "..", dir))
            append!(files, [joinpath(root, name) for name in sort(names) if endswith(name, ".jl")])
        end
    end
    return bytes2hex(sha256(string(VERSION, '\0', join(read.(files, String), '\0'))))
end

function check_prebuilt_examples(examples)
    for (category, exname) in examples
        base = joinpath(@__DIR__, "src", "examples", category, exname)
        marker = "$base.sha256"
        isfile("$base.md") && isfile(marker) &&
            read(marker, String) == example_fingerprint(category, exname) ||
            error("Missing or stale prebuilt example $category/$exname. Run docs/make.jl --examples-only first.")
    end
    return nothing
end

function prebuild_examples(examples; workers = 2)
    workers >= 1 || throw(ArgumentError("example_workers must be positive"))
    isempty(examples) && return nothing
    # Assign complete categories to a task to share compilation between related
    # examples. Tasks launch processes: cd, output capture and GLMakie are global
    # process state and cannot safely execute on concurrent Julia threads.
    batches = [Tuple{String, String}[] for _ in 1:min(workers, length(examples))]
    categories = unique(first.(examples))
    for (i, category) in enumerate(categories)
        append!(batches[mod1(i, length(batches))], filter(ex -> first(ex) == category, examples))
    end
    project = dirname(Base.active_project())
    worker = joinpath(@__DIR__, "example_worker.jl")
    logs = joinpath(@__DIR__, "example_logs")
    mkpath(logs)
    @sync for (i, batch) in enumerate(batches)
        isempty(batch) && continue
        @async begin
            args = ["$category/$exname" for (category, exname) in batch]
            cmd = `$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$project $worker $args`
            log = joinpath(logs, "worker-$i.log")
            @info "Starting example worker" i examples = args log
            open(log, "w") do io
                output = Pipe()
                try
                    proc = run(pipeline(addenv(cmd, "OPENBLAS_NUM_THREADS" => "1"),
                        stdout = output, stderr = output); wait = false)
                    close(output.in)
                    # Keep a complete worker log and show progress immediately in
                    # CI, including which example is running when a build stops.
                    for line in eachline(output)
                        println(io, line)
                        flush(io)
                        println("[example worker $i] ", line)
                        flush(stdout)
                    end
                    wait(proc)
                    success(proc) || error("Example worker $i failed with exit code $(proc.exitcode). See $log.")
                finally
                    close(output)
                end
            end
        end
    end
    check_prebuilt_examples(examples)
    return nothing
end

function replace_tags(content, subdir, exname)
    content_lines = split(content, "\n")
    for (i, line) in enumerate(content_lines)
        t = parse_tags(line)
        if !isnothing(t)
            content_lines[i] = tag_str(t)
            break
        end
    end
    content = join(content_lines, "\n")
    return content
end
