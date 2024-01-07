using JutulDarcy
using Jutul
using Literate
using Documenter

using DocumenterCitations
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

function build_jutul_darcy_docs(build_format = nothing; build_examples = true)
    DocMeta.setdocmeta!(JutulDarcy, :DocTestSetup, :(using JutulDarcy; using Jutul); recursive=true)
    DocMeta.setdocmeta!(Jutul, :DocTestSetup, :(using Jutul); recursive=true)

    ## Literate pass
    # Base directory
    jutul_dir = joinpath(dirname(pathof(JutulDarcy)), "..")
    # Convert examples as .jl files to markdown
    examples = [
        "Gravity segregation" => "two_phase_gravity_segregation",
        "Two-phase Buckley-Leverett" => "two_phase_buckley_leverett",
        "Gravity circulation with CPR preconditioner" => "two_phase_unstable_gravity",
        "Intro to wells" => "wells_intro",
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
        "[the JutulDarcy.jl GitHub repository](https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/$pth.jl)."
    end
    if build_examples
        for (ex, pth) in examples
            in_pth = joinpath(jutul_dir, "examples", "$pth.jl")
            out_dir = joinpath(jutul_dir, "docs", "src", "examples")
            push!(examples_markdown, ex => joinpath("examples", "$pth.md"))
            upd(content) = update_footer(content, pth)
            Literate.markdown(in_pth, out_dir, preprocess = upd)
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
            "Home" => "index.md",
            "Manual" => [
                "Getting started" =>"man/intro.md",
                "High-level functions" =>"man/highlevel.md",
                "Initial state" =>"man/state0.md",
                "Driving forces" =>"man/forces.md",
                "Supported physical systems" =>"man/systems.md",
                "Wells and controls" =>"man/wells.md",
                "Solving the equations" => "man/solution.md",
                "Primary variables" => "man/primary.md",
                "Secondary variables (properties)" => "man/secondary.md",
                "Parameters" => "man/parameters.md",
                "Input files" =>"man/input_files.md",
                "Visualization" =>"man/plotting.md",
                ],
            "Advanced usage" => [
                "Parallel solves with MPI" =>"man/mpi.md",
                "JutulDarcy.jl compiled app" =>"man/compiled.md"
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
        repo="github.com/sintefmath/JutulDarcy.jl",
        devbranch="main",
    )
end
##
# build_jutul_darcy_docs()

# ```@autodocs
# Modules = [JutulDarcy]
# ```
