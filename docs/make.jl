using JutulDarcy
using Jutul
using Literate
using Documenter
using OrderedCollections
using SHA

using DocumenterCitations
using DocumenterVitepress
using GLMakie
##
cd(@__DIR__)

# ENV["JUTULDARCY_RUN_VITEPRESS"] = 1
# ENV["JUTULDARCY_DOCS_EXAMPLES_SKIP"] = 1
# examples_to_build = ["validation_spe1"]

include("literate_examples.jl")

function build_jutul_darcy_docs(
        build_format = nothing;
        build_examples = true,
        build_docs = true,
        build_validation_examples = build_examples,
        build_notebooks = true,
        prebuilt_examples = false,
        example_workers = parse(Int, get(ENV, "JUTULDARCY_DOCS_EXAMPLE_WORKERS", "2")),
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
    if clean && !prebuilt_examples
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
    selected_examples = Tuple{String, String}[]
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
            push!(ex_dest, joinpath("examples", category, "$exname.md"))
            push!(selected_examples, (category, exname))
        end
    end
    write_tags()
    if prebuilt_examples
        check_prebuilt_examples(selected_examples)
    elseif get(ENV, "JUTULDARCY_DOCS_DRAFT_MODE", "0") == "1"
        # Preserve Documenter draft mode: convert pages without executing them.
        for (category, exname) in selected_examples
            build_literate_example(category, exname; execute = false)
        end
    else
        prebuild_examples(selected_examples; workers = example_workers)
    end
    examples_markdown = Any["examples/overview/example_overview.md"]
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
                    "man/intro.md",
                    "man/first_ex.md",
                    "man/python.md",
                    "FAQ" => "extras/faq.md",
                ],
                "Fundamentals" => [
                    "man/highlevel.md",
                    "man/basics/units.md",
                    "man/basics/input_files.md",
                    "man/basics/systems.md",
                    "man/basics/solution.md",
                ],
                "Detailed API" => [
                    "man/basics/wells.md",
                    "man/basics/adjoints.md",
                    "man/basics/primary.md",
                    "man/basics/secondary.md",
                    "man/basics/parameters.md",
                    "man/basics/plotting.md",
                    "man/basics/utilities.md",
                ],
                "Parallelism and compilation" => [
                    "man/advanced/parallel.md",
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
            "Models" => sort(validation_markdown),
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
            checkdocs_ignored_modules = [Jutul],
            draft = get(ENV, "JUTULDARCY_DOCS_DRAFT_MODE", "0") == "1"
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
# To use the standard documenter draft mode you can set
# ENV["JUTULDARCY_DOCS_DRAFT_MODE"] = 1
# This skips all examples, including the inline ones.
function main(args = ARGS)
    all(arg -> arg in ("--examples-only", "--docs-only"), args) || error("Unknown docs build argument: $args")
    length(args) <= 1 || error("Choose either --examples-only or --docs-only")
    if args == ["--examples-only"]
        build_jutul_darcy_docs(build_docs = false, build_notebooks = false, deploy = false)
        return
    end
    prebuilt = args == ["--docs-only"]
    if get(ENV, "JUTULDARCY_DOCS_EXAMPLES_SKIP", "0") == "1"
        # You can add a list of examples to build by defining examples_to_build.
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
        build_jutul_darcy_docs(prebuilt_examples = prebuilt)
    end
    if get(ENV, "JUTULDARCY_RUN_VITEPRESS", "0") == "1" && !Sys.iswindows()
        DocumenterVitepress.dev_docs("build")
    end
end

if get(ENV, "JUTULDARCY_DOCS_AUTOBUILD", "1") == "1"
    main()
end
