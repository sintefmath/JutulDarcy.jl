using JutulDarcy
using Jutul
using Literate
using Documenter
using OrderedCollections

using DocumenterCitations
using DocumenterVitepress
##
cd(@__DIR__)
function build_jutul_darcy_docs(
        build_format = nothing;
        build_examples = true,
        build_validation_examples = build_examples,
        build_notebooks = true,
        clean = true,
        deploy = true,
        use_vitepress = !Sys.iswindows()
    )
    DocMeta.setdocmeta!(JutulDarcy, :DocTestSetup, :(using JutulDarcy); recursive=true)
    DocMeta.setdocmeta!(Jutul, :DocTestSetup, :(using Jutul); recursive=true)
    bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

    ## Literate pass
    # Base directory
    jutul_dir = realpath(joinpath(@__DIR__, ".."))
    # Convert examples as .jl files to markdown
    examples = [
        # "Intro" => "intro_example",
        "Intro" => "two_phase_gravity_segregation",
        "Intro" => "two_phase_buckley_leverett",
        "Intro" => "wells_intro",
        "Intro" => "data_input_file",
        "Intro" => "intro_sensitivities",
        "General" => "five_spot_ensemble",
        "General" => "two_phase_unstable_gravity",
        "Workflow" => "co2_sloped",
        "Workflow" => "adding_new_wells",
        "Workflow" => "model_coarsening",
        "Workflow" => "hybrid_simulation_relperm",
        "Workflow" => "tracers_two_wells",
        "Data assimilation" => "cgnet_egg",
        "Data assimilation" => "optimize_simple_bl",
        "Compositional" => "co2_brine_2d_vertical",
        "Compositional" => "compositional_5components",
        "Discretizations" => "consistent_avgmpfa",
        "Discretizations" => "mpfa_weno_discretizations",
        "Properties" => "co2_props",
        "Properties" => "relperms",
        "Validation" => "validation_spe1",
        "Validation" => "validation_spe9",
        "Validation" => "validation_compositional",
        "Validation" => "validation_egg",
        "Validation" => "validation_olympus_1",
        "Validation" => "validation_norne_nohyst",
        "Validation" => "validation_mrst"
    ]
    validation_markdown = []
    examples_by_name = OrderedDict{String, Any}()
    function update_footer(content, pth)
        return content*"\n\n # ## Example on GitHub\n "*
        "# If you would like to run this example yourself, it can be downloaded from "*
        "the JutulDarcy.jl GitHub repository [as a script](https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/$pth.jl), "*
        "or as a [Jupyter Notebook](https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/$pth.ipynb)"
    end
    if clean
        for (ex, pth) in examples
            delpath = joinpath(@__DIR__, "src", "examples", "$pth.md")
            if isfile(delpath)
                println("Deleting generated example \"$ex\":\n\t$delpath")
                rm(delpath)
            else
                println("Did not find generated example \"$ex\", skipping removal:\n\t$delpath")
            end
        end
    end
    example_path(pth) = joinpath(jutul_dir, "examples", "$pth.jl")
    out_dir = joinpath(@__DIR__, "src", "examples")
    notebook_dir = joinpath(@__DIR__, "assets")
    for (category, pth) in examples
        in_pth = example_path(pth)
        is_validation = category == "Validation"
        if is_validation
            ex_dest = validation_markdown
            do_build = build_validation_examples
        else
            if !haskey(examples_by_name, category)
                examples_by_name[category] = []
            end
            ex_dest = examples_by_name[category]
            do_build = build_examples
        end
        if do_build
            push!(ex_dest, joinpath("examples", "$pth.md"))
            upd(content) = update_footer(content, pth)
            Literate.markdown(in_pth, out_dir, preprocess = upd)
        end
    end
    examples_markdown = []
    for (k, v) in pairs(examples_by_name)
        push!(examples_markdown, k => v)
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
        "Validation" => validation_markdown
    ]
    # for (k, subpages) in build_pages
    #     println("$k")
    #     @info "$k" subpages
    # end
    makedocs(;
        modules = [JutulDarcy, Jutul],
        authors = "Olav Møyner <olav.moyner@sintef.no> and contributors",
        repo = "https://github.com/sintefmath/JutulDarcy.jl/blob/{commit}{path}#{line}",
        sitename = "Documentation | JutulDarcy.jl",
        warnonly = false,
        checkdocs = :exports,
        plugins = [bib],
        format = build_format,
        pages = build_pages,
    )
    if build_notebooks
        # Subfolder of final site build folder
        notebook_dir = joinpath(@__DIR__, "build", "final_site", "notebooks")
        mkpath(notebook_dir)
        for (ex, pth) in examples
            in_pth = example_path(pth)
            @info "$ex Writing notebook to $notebook_dir"
            Literate.notebook(in_pth, notebook_dir, execute = false)
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
##
# build_jutul_darcy_docs(
#     build_examples = false,
#     build_validation_examples = false,
#     build_notebooks = false,
#     deploy = false
# )
build_jutul_darcy_docs()

# ```@autodocs
# Modules = [JutulDarcy]
# ```
