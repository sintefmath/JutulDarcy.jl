using JutulDarcy
using Jutul
using Literate
using Documenter

DocMeta.setdocmeta!(JutulDarcy, :DocTestSetup, :(using JutulDarcy; using Jutul); recursive=true)
## Literate pass
# Base directory
jutul_dir = joinpath(dirname(pathof(JutulDarcy)), "..")
# Convert examples as .jl files to markdown
examples = ["Intro to wells" => "wells_intro"]
examples_markdown = ["Getting started" => "examples/intro.md"]
for (ex, pth) in examples
    in_pth = joinpath(jutul_dir, "examples", "$pth.jl")
    out_dir = joinpath(jutul_dir, "docs", "src", "examples")
    push!(examples_markdown, ex => joinpath("examples", "$pth.md"))
    Literate.markdown(in_pth, out_dir)
end
## Docs
makedocs(;
    modules=[JutulDarcy, Jutul],
    authors="Olav Møyner <olav.moyner@sintef.no> and contributors",
    repo="https://github.com/sintefmath/JutulDarcy.jl/blob/{commit}{path}#{line}",
    sitename="JutulDarcy.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://sintefmath.github.io/JutulDarcy.jl",
        edit_link="main",
        assets=String[],
    ),
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
