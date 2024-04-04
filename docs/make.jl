using JutulDarcy
using Jutul
using Literate
using Documenter

using DocumenterCitations
##
function build_jutul_darcy_docs(build_format = nothing; build_examples = true, build_notebooks = build_examples, clean = true)
    DocMeta.setdocmeta!(JutulDarcy, :DocTestSetup, :(using JutulDarcy; using Jutul); recursive=true)
    DocMeta.setdocmeta!(Jutul, :DocTestSetup, :(using Jutul); recursive=true)
    bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

    ## Literate pass
    # Base directory
    jutul_dir = realpath(joinpath(@__DIR__, ".."))
    # Convert examples as .jl files to markdown
    examples = [
        "Intro" => "intro",
        "Gravity segregation" => "two_phase_gravity_segregation",
        "Two-phase Buckley-Leverett" => "two_phase_buckley_leverett",
        "Gravity circulation with CPR preconditioner" => "two_phase_unstable_gravity",
        "Intro to wells" => "wells_intro",
        "CO2 injection in saline aquifer" => "co2_sloped",
        "Quarter-five-spot with variation" => "five_spot_ensemble",
        "Intro to compositional flow" => "co2_brine_2d_vertical",
        "Compositional with five components" => "compositional_5components",
        "Parameter optimization of Buckley-Leverett" => "optimize_simple_bl",
        "Validation of reservoir simulator" => "mrst_validation"
    ]
    examples_markdown = ["Getting started" => "examples/intro.md"]
    function update_footer(content, pth)
        return content*"\n\n # ## Example on GitHub\n "*
        "# If you would like to run this example yourself, it can be downloaded from "*
        "the JutulDarcy.jl GitHub repository [as a script](https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/$pth.jl), "*
        "or as a [Notebook](https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/docs/notebooks/$pth.ipynb)"
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
    for (ex, pth) in examples
        in_pth = example_path(pth)
        if build_examples
            out_dir = joinpath(@__DIR__, "src", "examples")
            push!(examples_markdown, ex => joinpath("examples", "$pth.md"))
            upd(content) = update_footer(content, pth)
            Literate.markdown(in_pth, out_dir, preprocess = upd)
        end
        if build_notebooks
            out_dir_notebooks = joinpath(@__DIR__, "notebooks")
            Literate.notebook(in_pth, out_dir_notebooks)
        end
    end
    ## Docs
    if isnothing(build_format)
        build_format = Documenter.HTML(;
            prettyurls=get(ENV, "CI", "false") == "true",
            canonical="https://sintefmath.github.io/JutulDarcy.jl",
            edit_link="main",
            assets=String["assets/citations.css"],
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
                "man/highlevel.md",
                ],
            "Manual" => [
                "man/basics/input_files.md",
                "man/basics/forces.md",
                "man/basics/systems.md",
                "man/basics/wells.md",
                "man/basics/solution.md",
                "man/basics/primary.md",
                "man/basics/secondary.md",
                "man/basics/parameters.md",
                "man/basics/plotting.md",
                ],
            "Advanced usage" => [
                "man/advanced/mpi.md",
                "man/advanced/compiled.md"
            ],
            "Examples" => examples_markdown,
            "Reference" => [
                # "Internals" => "ref/internals.md",
                "Jutul functions" => "ref/jutul.md"
            ],
            "Additional information "=> [
                "References" => "extras/refs.md",
                "FAQ" => "extras/faq.md"
            ]
        ],
    )

    deploydocs(;
        repo="github.com/sintefmath/JutulDarcy.jl.git",
        devbranch="main",
    )
end
##
build_jutul_darcy_docs()

# ```@autodocs
# Modules = [JutulDarcy]
# ```
