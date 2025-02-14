using JutulDarcy
using Jutul
using Literate
using Documenter
using OrderedCollections

using DocumenterCitations
using DocumenterVitepress
##
cd(@__DIR__)

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

function update_footer(content, subdir, exname)
    return content*"\n\n # ## Example on GitHub\n "*
    "# If you would like to run this example yourself, it can be downloaded from "*
    "the JutulDarcy.jl GitHub repository [as a script](https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/$subdir/$exname.jl), "*
    "or as a [Jupyter Notebook](https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/$subdir/$exname.ipynb)"
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
    # Base directory
    jutul_dir = realpath(joinpath(@__DIR__, ".."))
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
    example_path(cname, pth) = joinpath(jutul_dir, "examples", cname, "$pth.jl")
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
            in_pth = example_path(category, exname)
            push!(ex_dest, joinpath("examples", category, "$exname.md"))
            upd(content) = update_footer(content, category, exname)
            Literate.markdown(in_pth, joinpath(out_dir, category), preprocess = upd)
        end
    end
    examples_markdown = []
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
                in_pth = example_path(category, exname)
                ex_notebook_dir = joinpath(notebook_dir, category)
                @info "$exname Writing notebook to $ex_notebook_dir"
                Literate.notebook(in_pth, ex_notebook_dir, execute = false)
            end
        end
    end
    if deploy
        deploydocs(;
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
