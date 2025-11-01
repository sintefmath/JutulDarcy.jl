#!/usr/bin/env julia

# Simple script to generate an example overview for the PDF
# This doesn't require any dependencies

function get_example_files()
    basepth = joinpath(@__DIR__, "..", "examples")
    examples = Dict{String, Vector{String}}()
    
    categories = ["introduction", "workflow", "data_assimilation", "geothermal", 
                  "compositional", "discretization", "properties", "validation"]
    
    for category in categories
        catpath = joinpath(basepth, category)
        if isdir(catpath)
            examples[category] = String[]
            for file in readdir(catpath)
                if endswith(file, ".jl")
                    filename = first(splitext(file))
                    push!(examples[category], filename)
                end
            end
            sort!(examples[category])
        end
    end
    
    return examples, categories
end

function category_title(cat)
    return titlecase(replace(cat, "_" => " "))
end

function write_example_overview()
    examples, categories = get_example_files()
    
    outdir = joinpath(@__DIR__, "src", "examples", "overview")
    mkpath(outdir)
    outpath = joinpath(outdir, "example_overview.md")
    
    open(outpath, "w") do io
        println(io, "# Example Overview\n")
        println(io, "JutulDarcy.jl comes with a number of examples that illustrate different features of the simulator.\n")
        println(io, "The examples are organized by category. For the full interactive examples with code and outputs, please visit the online documentation at https://sintefmath.github.io/JutulDarcy.jl/\n")
        
        for category in categories
            if haskey(examples, category) && !isempty(examples[category])
                println(io, "## $(category_title(category))\n")
                
                # Add category descriptions
                if category == "introduction"
                    println(io, "Basic examples that illustrate fundamental features of JutulDarcy.jl.\n")
                elseif category == "workflow"
                    println(io, "Examples demonstrating complete workflows and advanced use cases.\n")
                elseif category == "data_assimilation"
                    println(io, "Examples of history matching, optimization, and sensitivity analysis.\n")
                elseif category == "geothermal"
                    println(io, "Geothermal reservoir simulation examples.\n")
                elseif category == "compositional"
                    println(io, "Compositional flow and multi-component examples.\n")
                elseif category == "discretization"
                    println(io, "Examples showing different discretization schemes.\n")
                elseif category == "properties"
                    println(io, "Examples focusing on fluid properties and relationships.\n")
                elseif category == "validation"
                    println(io, "Validation cases comparing with other simulators and benchmarks.\n")
                end
                
                for ex in examples[category]
                    # Convert underscores to readable format
                    title = titlecase(replace(ex, "_" => " "))
                    println(io, "- **$title** (`$ex.jl`)")
                end
                println(io, "")
            end
        end
        
        println(io, "\n---\n")
        println(io, "*Note: The full examples with code, plots, and detailed explanations are available in the [online documentation](https://sintefmath.github.io/JutulDarcy.jl/dev/). ")
        println(io, "You can also find the example scripts in the `examples/` directory of the repository.*")
    end
    
    println("Generated example overview at: $outpath")
end

# Run the function
write_example_overview()
