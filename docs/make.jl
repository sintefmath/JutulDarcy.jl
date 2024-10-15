using JutulDarcy
using Jutul
using Literate
using Documenter

using DocumenterCitations
using DocumenterVitepress
##
cd(@__DIR__)
function build_jutul_darcy_docs(build_format = nothing;
        build_examples = true,
        build_validation_examples = build_examples,
        build_notebooks = true,
        clean = true,
        deploy = true
    )
    DocMeta.setdocmeta!(JutulDarcy, :DocTestSetup, :(using JutulDarcy; using Jutul); recursive=true)
    DocMeta.setdocmeta!(Jutul, :DocTestSetup, :(using Jutul); recursive=true)
    bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

    ## Literate pass
    # Base directory
    jutul_dir = realpath(joinpath(@__DIR__, ".."))
    # Convert examples as .jl files to markdown
    examples = [
        # "Intro" => "intro_example",
        "Intro: Gravity segregation" => "two_phase_gravity_segregation",
        "Intro: Two-phase Buckley-Leverett" => "two_phase_buckley_leverett",
        "Intro: Wells" => "wells_intro",
        "Intro: Simulating Eclipse/DATA input files" => "data_input_file",
        "Intro: Sensitivities in JutulDarcy" => "intro_sensitivities",
        "Intro: Compositional flow" => "co2_brine_2d_vertical",
        "Quarter-five-spot with variation" => "five_spot_ensemble",
        "Gravity circulation with CPR preconditioner" => "two_phase_unstable_gravity",
        "CO2 property calculations" => "co2_props",
        "Relative permeability" => "relperms",
        "Model coarsening" => "model_coarsening",
        "History matching a coarse model - CGNet" => "cgnet_egg",
        "CO2 injection in saline aquifer" => "co2_sloped",
        "Compositional with five components" => "compositional_5components",
        "Parameter matching of Buckley-Leverett" => "optimize_simple_bl",
        "Validation: SPE1" => "validation_spe1",
        "Validation: SPE9" => "validation_spe9",
        "Validation: Compositional" => "validation_compositional",
        "Validation: Egg" => "validation_egg",
        "Validation: OLYMPUS 1" => "validation_olympus_1",
        "Validation: Norne" => "validation_norne_nohyst",
        "Validation: MRST input files" => "validation_mrst"
    ]
    examples_markdown = []
    validation_markdown = []
    intros_markdown = []
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
    for (ex, pth) in examples
        in_pth = example_path(pth)
        is_validation = startswith(ex, "Validation:")
        is_intro = startswith(ex, "Intro: ")
        is_example = !(is_intro || is_validation)
        if is_validation
            ex_dest = validation_markdown
            do_build = build_validation_examples
        else
            if is_intro
                ex_dest = intros_markdown
            else
                ex_dest = examples_markdown
            end
            do_build = build_examples
        end
        if do_build
            push!(ex_dest, ex => joinpath("examples", "$pth.md"))
            upd(content) = update_footer(content, pth)
            Literate.markdown(in_pth, out_dir, preprocess = upd)
        end
    end
    ## Docs
    if isnothing(build_format)
        # Old Documenter code in case we want to go back.
        # build_format = Documenter.HTML(;
        #     prettyurls=get(ENV, "CI", "false") == "true",
        #     canonical="https://sintefmath.github.io/JutulDarcy.jl",
        #     edit_link="main",
        #     size_threshold_ignore = ["ref/jutul.md", "docstrings.md"],
        #     assets=String["assets/citations.css"],
        # )
        build_format = DocumenterVitepress.MarkdownVitepress(
            repo = "https://github.com/sintefmath/JutulDarcy.jl",
        )
    end
    makedocs(;
        modules=[JutulDarcy, Jutul],
        authors="Olav Møyner <olav.moyner@sintef.no> and contributors",
        repo="https://github.com/sintefmath/JutulDarcy.jl/blob/{commit}{path}#{line}",
        sitename="JutulDarcy.jl",
        warnonly = true,
        plugins=[bib],
        format=build_format,
        pages=[
            "Introduction" => [
                "JutulDarcy.jl" => "index.md",
                "Getting started" =>"man/intro.md",
                "Your first JutulDarcy.jl simulation" => "man/first_ex.md",
                "FAQ" => "extras/faq.md"
            ],
            "Manual" => [
                    "man/highlevel.md",
                    "man/basics/input_files.md",
                    "man/basics/forces.md",
                    "man/basics/systems.md",
                    "man/basics/wells.md",
                    "man/basics/solution.md",
                    "man/basics/primary.md",
                    "man/basics/secondary.md",
                    "man/basics/parameters.md",
                    "man/basics/plotting.md",
                    "man/basics/utilities.md",
                    ],
            "Further reading" => [
                "man/advanced/mpi.md",
                "man/advanced/compiled.md",
                "Jutul functions" => "ref/jutul.md",
                "Bibliography" => "extras/refs.md"
            ],
            "Examples: Introduction" => intros_markdown,
            "Examples: Usage" => examples_markdown,
            "Examples: Validation" => validation_markdown
        ],
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
##
# build_jutul_darcy_docs(build_examples = false, build_validation_examples = false)
build_jutul_darcy_docs()

# ```@autodocs
# Modules = [JutulDarcy]
# ```
