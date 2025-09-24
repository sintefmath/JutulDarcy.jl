using JutulDarcy
using Jutul
using Literate
using Documenter
using OrderedCollections

using DocumenterCitations
using DocumenterVitepress
using GLMakie
##
cd(@__DIR__)

# ENV["JUTULDARCY_RUN_VITEPRESS"] = 1
# ENV["JUTULDARCY_DOCS_EXAMPLES_SKIP"] = 1
# examples_to_build = ["validation_spe1"]

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

example_tags(example_path_jl("validation", "validation_spe1"))

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
    descr["Geothermal"] = "Examples that simulate recovery and/or storage of heat in the subsurface. See also the dedicated [Fimbul.jl](https://sintefmath.github.io/Fimbul.jl/dev/) module for geothermal simulation with JutulDarcy.jl."

    out = OrderedDict()
    colors = to_colormap(:tab20)
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
    return out
end

function write_tags()
    tags = all_tags()
    outpth = joinpath(@__DIR__, "src", "example_overview.md")
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
                    exlink = joinpath("examples", category, "$exname.md")
                    println(io, "1. [$exname]($exlink) (in $category)")
                end
                println(io, "\n")
            end
        end
    end
    println("Wrote tags to $outpth")
end

write_tags()

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
    "the JutulDarcy.jl GitHub repository [as a script](https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/$subdir/$exname.jl), "*
    "or as a [Jupyter Notebook](https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/$subdir/$exname.ipynb)"
end

function update_footer(content, subdir, exname)
    info_footer = example_info_footer(subdir, exname)
    gc_footer = post_run_variables_gc()
    start, stop = timer_str()
    new_content = string(start, content, info_footer, stop, gc_footer)
    # print(new_content)
    return new_content
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

function build_jutul_darcy_docs(
        build_format = nothing;
        build_examples = true,
        build_docs = true,
        build_validation_examples = build_examples,
        build_notebooks = true,
        examples_explicit_list = missing,
        skip_examples = String[],
        clean = true,
        deploy = true,
        use_vitepress = !Sys.iswindows()
    )
    if examples_explicit_list isa String
        examples_explicit_list = [examples_explicit_list]
    end
    examples_explicit_list::Union{Vector{String}, Missing}
    has_explicit_list = !ismissing(examples_explicit_list)
    if has_explicit_list
        @info "Building only examples as examples_explicit_list was specified" examples_explicit_list
    end
    DocMeta.setdocmeta!(JutulDarcy, :DocTestSetup, :(using JutulDarcy); recursive=true)
    DocMeta.setdocmeta!(Jutul, :DocTestSetup, :(using Jutul); recursive=true)
    bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

    ## Literate pass
    # Convert examples as .jl files to markdown
    examples = get_example_paths(check_empty = !has_explicit_list)
    validation_markdown = []
    examples_by_name = OrderedDict{String, Any}()
    if clean
        for (category, example_set) in pairs(examples)
            for ex in example_set
                delpath = joinpath(@__DIR__, "src", "examples", category, "$ex.md")
                if isfile(delpath)
                    println("Deleting generated example \"$ex\":\n\t$delpath")
                    rm(delpath)
                else
                    println("Did not find generated example \"$ex\", skipping removal:\n\t$delpath")
                end
            end
        end
    end
    out_dir = joinpath(@__DIR__, "src", "examples")
    notebook_dir = joinpath(@__DIR__, "assets")
    for (category, example_set) in pairs(examples)
        if category == "validation"
            ex_dest = validation_markdown
            do_build = build_validation_examples
        else
            ex_dest = []
            examples_by_name[category] = ex_dest
            do_build = build_examples
        end
        for exname in example_set
            if has_explicit_list
                if exname in examples_explicit_list
                    jutul_message("Examples", "$category/$exname added to build from explicit list.", color = :green)
                else
                    jutul_message("Examples", "$category/$exname not in explicit list, skipping.", color = :yellow)
                    continue
                end
            elseif do_build && !(exname in skip_examples)
                jutul_message("Examples", "$category/$exname was added.", color = :green)
            else
                jutul_message("Examples", "$category/$exname was skipped.", color = :blue)
                continue
            end
            in_pth = example_path_jl(category, exname)
            push!(ex_dest, joinpath("examples", category, "$exname.md"))
            upd(content) = update_footer(content, category, exname)
            fixt(content) = replace_tags(content, category, exname)
            Literate.markdown(in_pth, joinpath(out_dir, category), preprocess = upd, postprocess = fixt)
        end
    end
    examples_markdown = Any["example_overview.md"]
    for (k, v) in pairs(examples_by_name)
        push!(examples_markdown, dir_to_doc_name(k) => v)
    end

    ## Docs
    if isnothing(build_format)
        if use_vitepress
            build_format = DocumenterVitepress.MarkdownVitepress(
                repo = "https://github.com/sintefmath/JutulDarcy.jl",
            )
        else
            build_format = Documenter.HTML(;
                prettyurls=get(ENV, "CI", "false") == "true",
                canonical="https://sintefmath.github.io/JutulDarcy.jl",
                edit_link="main",
                size_threshold_ignore = [
                    "ref/jutul.md",
                    "docstrings.md",
                    "man/first_ex.md"
                ],
                assets=String["assets/citations.css"],
            )
        end
    end
    build_pages = [
        "Manual" => [
                "Introduction" => [
                    "JutulDarcy.jl" => "index.md",
                    "Getting started" =>"man/intro.md",
                    "Your first JutulDarcy.jl simulation" => "man/first_ex.md",
                    "FAQ" => "extras/faq.md",
                ],
                "Fundamentals" => [
                    "man/highlevel.md",
                    "man/basics/input_files.md",
                    "man/basics/systems.md",
                    "man/basics/solution.md",
                ],
                "Detailed API" => [
                    "man/basics/forces.md",
                    "man/basics/wells.md",
                    "man/basics/primary.md",
                    "man/basics/secondary.md",
                    "man/basics/parameters.md",
                    "man/basics/plotting.md",
                    "man/basics/utilities.md",
                ],
                "Parallelism and compilation" => [
                    "man/advanced/mpi.md",
                    "man/advanced/gpu.md",
                    "man/advanced/compiled.md"
                ],
                "References" => [
                    "man/basics/package.md",
                    "extras/paper_list.md",
                    "Jutul functions" => "ref/jutul.md",
                    "Bibliography" => "extras/refs.md"
                ],
            ],
        "Examples" => examples_markdown,
        "Validation" => [
            "man/validation.md",
            "Models" => validation_markdown,
        ]
    ]
    # for (k, subpages) in build_pages
    #     println("$k")
    #     @info "$k" subpages
    # end
    if build_docs
        makedocs(;
            modules = [JutulDarcy, Jutul],
            authors = "Olav Møyner <olav.moyner@sintef.no> and contributors",
            repo = "https://github.com/sintefmath/JutulDarcy.jl/blob/{commit}{path}#{line}",
            warnonly = false,
            sitename = "JutulDarcy.jl",
            checkdocs = :exports,
            plugins = [bib],
            format = build_format,
            pages = build_pages,
        )
    end
    if build_notebooks
        # Subfolder of final site build folder
        notebook_dir = joinpath(@__DIR__, "build", "final_site", "notebooks")
        mkpath(notebook_dir)
        for (category, example_set) in pairs(examples)
            for exname in example_set
                in_pth = example_path_jl(category, exname)
                ex_notebook_dir = joinpath(notebook_dir, category)
                @info "$exname Writing notebook to $ex_notebook_dir"
                Literate.notebook(in_pth, ex_notebook_dir, execute = false)
            end
        end
    end
    if deploy
        DocumenterVitepress.deploydocs(;
            repo="github.com/sintefmath/JutulDarcy.jl.git",
            devbranch="main",
            target = "build", # this is where Vitepress stores its output
            branch = "gh-pages",
            push_preview = true
        )
    end
end
# To preview, go to the docs folder and run:
# # DocumenterVitepress.dev_docs("build")
# To only build some examples you can set
# ENV["JUTULDARCY_DOCS_EXAMPLES_SKIP"] = 1
# You can also enable build after (Linux only):
# ENV["JUTULDARCY_RUN_VITEPRESS"] = 1
if get(ENV, "JUTULDARCY_DOCS_EXAMPLES_SKIP", "0") == "1"
    # You can add a list of examples to build by running
    # examples_to_build = ["geothermal_1well"]
    if isdefined(Main, :examples_to_build)
        ex_list = examples_to_build
    else
        ex_list = missing
    end
    build_jutul_darcy_docs(
        build_examples = false,
        build_validation_examples = false,
        build_notebooks = false,
        examples_explicit_list = ex_list,
        deploy = false
    )
else
    build_jutul_darcy_docs()
end

if get(ENV, "JUTULDARCY_RUN_VITEPRESS", "0") == "1" && !Sys.iswindows()
    DocumenterVitepress.dev_docs("build")
end
