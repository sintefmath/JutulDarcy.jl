using Test, Literate, OrderedCollections, SHA, Documenter
include("../literate_examples.jl")

@testset "Prebuilt Literate examples" begin
    repo = dirname(dirname(dirname(@__FILE__)))
    inputs = [mktempdir(joinpath(repo, "examples"); prefix = "docs_test_") for _ in 1:2]
    categories = basename.(inputs)
    outputs = [joinpath(repo, "docs", "src", "examples", category) for category in categories]
    examples = [(category, "index") for category in categories]
    try
        for (i, input) in enumerate(inputs)
            counter = repr(joinpath(input, "runs.txt"))
            source = """
            # # Example execution fixture
            # <tags: Properties>
            open($counter, "a") do io # hide
                println(io, "run") # hide
            end # hide
            hidden_value = 4 # hide
            #-
            hidden_value + 3
            #-
            display("explicit display output")
            #-
            println("captured stdout")
            #-
            struct FixtureImage end # hide
            Base.show(io::IO, ::MIME"image/svg+xml", ::FixtureImage) = print(io, "<svg xmlns='http://www.w3.org/2000/svg' width='10' height='10'></svg>") # hide
            FixtureImage()
            #-
            struct FixtureHTMLImage end # hide
            Base.show(io::IO, ::MIME"text/html", ::FixtureHTMLImage) = print(io, "<img src='placeholder'>") # hide
            FixtureHTMLImage()
            #-
            Dict(:large_table => collect(1:10000))
            # ## A second large result
            collect(1:10000)'
            """
            source = replace(source, "src='placeholder'" => "src=\\\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVQIHWP4z8DwHwAFgAI/ScLbtAAAAABJRU5ErkJggg==\\\"")
            write(joinpath(input, "index.jl"), i == 1 ? replace(source, "\n" => "\r\n") : source)
        end
        prebuild_examples(examples; workers = 2)
        for output in outputs
            md = read(joinpath(output, "index.md"), String)
            @test !occursin("@example", md)
            @test !occursin("hidden_value = 4", md)
            @test occursin("7", md)
            @test occursin("explicit display output", md)
            @test occursin("captured stdout", md)
            @test occursin("Example on GitHub", md)
            @test occursin("EditURL", md)
            @test occursin("<ExampleTag", md)
            @test any(endswith(".svg"), readdir(output))
            @test !occursin("data:image/png;base64", md)
            @test occursin("src=\"./index-inline-1.png\"", md)
            @test read(joinpath(output, "index-inline-1.png"))[1:8] == UInt8[137, 80, 78, 71, 13, 10, 26, 10]
            @test occursin("output truncated", md)
            @test length(findall("output truncated", md)) == 2
            @test occursin("Dict(:large_table => collect(1:10000))", md)
            @test sizeof(md) < 10_000
        end
        check_prebuilt_examples(examples)
        # Documenter consumes the evaluated Markdown without running the scripts again.
        mktempdir() do build
            makedocs(; root = joinpath(repo, "docs"), source = first(outputs), build,
                sitename = "Prebuild smoke test", remotes = nothing, doctest = false,
                checkdocs = :none, format = Documenter.HTML(prettyurls = false, edit_link = nothing))
            @test isfile(joinpath(build, "index.html"))
        end
        @test all(input -> readlines(joinpath(input, "runs.txt")) == ["run"], inputs)
        open(joinpath(first(inputs), "index.jl"), "a") do io
            println(io, "# changed source")
        end
        @test_throws ErrorException check_prebuilt_examples(examples)
        write(joinpath(first(inputs), "index.jl"), "error(\"example failure\")\n")
        @test_throws CompositeException prebuild_examples(examples[1:1]; workers = 1)
        @test !isfile(joinpath(first(outputs), "index.sha256"))
    finally
        for dir in vcat(inputs, outputs)
            @assert startswith(normpath(dir), joinpath(repo, "examples", "")) ||
                startswith(normpath(dir), joinpath(repo, "docs", "src", "examples", ""))
            rm(dir; recursive = true, force = true)
        end
    end
end
