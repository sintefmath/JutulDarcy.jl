using JutulDarcy
using Jutul
using Literate
using Documenter

function build_jutul_darcy_docs(build_format = nothing; build_examples = true)
    DocMeta.setdocmeta!(JutulDarcy, :DocTestSetup, :(using JutulDarcy; using Jutul); recursive=true)
    ## Literate pass
    # Base directory
    jutul_dir = joinpath(dirname(pathof(JutulDarcy)), "..")
    # Convert examples as .jl files to markdown
    examples = [
        "Two-phase Buckley-Leverett" => "two_phase_buckley_leverett",
        "Intro to wells" => "wells_intro"
    ]
    examples_markdown = ["Getting started" => "examples/intro.md"]
    if build_examples
        for (ex, pth) in examples
            in_pth = joinpath(jutul_dir, "examples", "$pth.jl")
            out_dir = joinpath(jutul_dir, "docs", "src", "examples")
            push!(examples_markdown, ex => joinpath("examples", "$pth.md"))
            Literate.markdown(in_pth, out_dir)
        end
    end
    ## Docs
    if isnothing(build_format)
        build_format = Documenter.HTML(;
            prettyurls=get(ENV, "CI", "false") == "true",
            canonical="https://sintefmath.github.io/JutulDarcy.jl",
            edit_link="main",
            assets=String[],
        )
    end
    makedocs(;
        modules=[JutulDarcy, Jutul],
        authors="Olav Møyner <olav.moyner@sintef.no> and contributors",
        repo="https://github.com/sintefmath/JutulDarcy.jl/blob/{commit}{path}#{line}",
        sitename="JutulDarcy.jl",
        format=build_format,
        pages=[
            "Home" => "index.md",
            "Examples" => examples_markdown,
            "Usage" => "usage.md",
            "Internals" => "internals.md"
        ],
    )

    deploydocs(;
        repo="github.com/sintefmath/JutulDarcy.jl",
        devbranch="main",
    )
end

build_jutul_darcy_docs()
