"""
Embedded fracture models on an unchanged three-dimensional matrix mesh.
"""
module EFM
    using Jutul, JutulDarcy, LinearAlgebra, StaticArrays
    import JutulDarcy: FractureMatrixConnection

    export AbstractEmbeddedFracture, PolygonFracture, RectangleFracture
    export EllipticalFracture, EllipsoidFracture, fracture_polygon
    export embedded_fracture_domain, embedded_fracture_mesh
    export combined_fracture_domain, setup_embedded_fracture_model
    export EmbeddedFractureGeometry, EmbeddedFractureMesh, EmbeddedReservoirMesh
    export discrete_fracture_mesh, discrete_fracture_domain
    export FractureCollection

    include("fractures.jl")
    include("collection.jl")
    include("geometry.jl")
    include("spatial.jl")
    include("mesh.jl")
    include("utils.jl")
    include("unstructured.jl")
    include("preprocessing.jl")
    include("domains.jl")
    include("discrete.jl")
    include("wells.jl")
    include("models.jl")
    include("unified.jl")
    include("unified_dfm.jl")
    include("show.jl")
end
