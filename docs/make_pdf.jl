using JutulDarcy
using Jutul
using Documenter
using DocumenterCitations
using DocumenterLaTeX
using OrderedCollections

##
cd(@__DIR__)

# Set up document metadata
DocMeta.setdocmeta!(JutulDarcy, :DocTestSetup, :(using JutulDarcy); recursive=true)
DocMeta.setdocmeta!(Jutul, :DocTestSetup, :(using Jutul); recursive=true)

# Set up bibliography
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

# Define the documentation structure
# Note: We're building a simplified structure for PDF that doesn't include
# dynamically generated examples from Literate.jl
build_pages = [
    "Manual" => [
        "Introduction" => [
            "JutulDarcy.jl" => "index_pdf.md",
            "Getting started" => "man/intro.md",
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
    "Validation" => [
        "man/validation.md",
    ]
]

# Configure LaTeX output
latex_format = DocumenterLaTeX.LaTeX(
    platform = "docker"  # Use docker platform for better compatibility
)

# Build the documentation
makedocs(;
    modules = [JutulDarcy, Jutul],
    authors = "Olav Møyner <olav.moyner@sintef.no> and contributors",
    repo = "https://github.com/sintefmath/JutulDarcy.jl/blob/{commit}{path}#{line}",
    sitename = "JutulDarcy.jl",
    format = latex_format,
    pages = build_pages,
    plugins = [bib],
    checkdocs = :exports,
    warnonly = true,  # Don't fail on warnings for PDF build
)

println("\n" * "="^80)
println("PDF documentation build complete!")
println("The LaTeX files are in: $(joinpath(@__DIR__, "build"))")
println("To compile the PDF, you'll need a LaTeX distribution installed.")
println("="^80)
