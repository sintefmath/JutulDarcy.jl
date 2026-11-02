using Literate, OrderedCollections, GLMakie, SHA
include("literate_examples.jl")
GLMakie.activate!(visible = false)

for example in ARGS
    category, exname = split(example, '/'; limit = 2)
    build_literate_example(String(category), String(exname))
    GLMakie.closeall()
    GC.gc()
end
