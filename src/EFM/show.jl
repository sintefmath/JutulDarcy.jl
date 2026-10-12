function show_fracture_count(io, count, label; plural = string(label, "s"))
    print(io, count, " ")
    if count == 1
        print(io, label)
    else
        print(io, plural)
    end
end

function Base.show(io::IO, fractures::FractureCollection)
    print(io, "FractureCollection(")
    show_fracture_count(io, length(fractures), "fracture")
    print(io, ", ")
    show_fracture_count(io, length(fractures.groups), "type group")
    print(io, ")")
end

function show_fracture_groups(io, fractures)
    for group in fractures.groups
        print(io, "\n  ", length(group), " × ")
        show(io, eltype(group))
    end
end

function Base.show(io::IO, ::MIME"text/plain", fractures::FractureCollection)
    show(io, fractures)
    show_fracture_groups(io, fractures)
end

function Base.show(io::IO, geometry::EmbeddedFractureGeometry{T}) where T
    print(io, "EmbeddedFractureGeometry{", T, "}(")
    show_fracture_count(io, length(geometry.patches), "patch")
    print(io, ", ")
    show_fracture_count(io, length(geometry.connection_types), "fracture connection")
    print(io, ", ")
    show_fracture_count(io, length(geometry.matrix_cells), "matrix connection")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", geometry::EmbeddedFractureGeometry)
    show(io, geometry)
    print(io, "\n  ")
    show_fracture_count(io, length(geometry.frames), "fracture plane")
    print(io, ", ")
    show_fracture_count(io, length(geometry.intersections), "intersection segment")
    print(io, "\n  Clipping tolerance: ", geometry.tolerance)
    print(io, "\n  ")
    show_fracture_count(io, count(grid -> !isnothing(grid), geometry.grids), "independent fracture grid")
end

function Base.show(io::IO, mesh::EmbeddedFractureMesh{T}) where T
    print(io, "EmbeddedFractureMesh{", T, "}(")
    show_fracture_count(io, length(mesh.fractures), "fracture")
    print(io, ", ")
    show_fracture_count(io, number_of_cells(mesh), "cell")
    print(io, ", ")
    show_fracture_count(io, number_of_faces(mesh), "connection")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", mesh::EmbeddedFractureMesh)
    show(io, mesh)
    print(io, "\n  Original matrix: ")
    show_fracture_count(io, number_of_cells(mesh.matrix_mesh), "cell")
    show_fracture_groups(io, mesh.fractures)
    print(io, "\n  ")
    show_fracture_count(io, length(mesh.geometry.matrix_cells), "matrix–fracture connection")
    print(io, "; ")
    show_fracture_count(io, length(mesh.geometry.intersections), "intersection segment")
end

function Base.show(io::IO, mesh::EmbeddedReservoirMesh{T}) where T
    print(io, "EmbeddedReservoirMesh{", T, "}(")
    show_fracture_count(io, number_of_cells(mesh.matrix_mesh), "matrix cell")
    print(io, ", ")
    show_fracture_count(io, number_of_cells(mesh.fracture_mesh), "fracture cell")
    print(io, ", ")
    show_fracture_count(io, number_of_faces(mesh), "connection")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", mesh::EmbeddedReservoirMesh)
    show(io, mesh)
    if mesh.fracture_mesh isa EmbeddedFractureMesh
        show_fracture_groups(io, mesh.fracture_mesh.fractures)
    else
        print(io, "\n  Cut-cell DFM fracture network")
    end
end

function Base.show(io::IO, fracture::RectangleFracture{T}) where T
    print(io, "RectangleFracture{", T, "}(center = ", Tuple(fracture.center),
        ", lengths = ", fracture.lengths, ", aperture = ", fracture.properties.aperture)
    show_fracture_mesh_options(io, fracture.mesh_options)
    print(io, ")")
end

function Base.show(io::IO, fracture::EllipticalFracture{T}) where T
    print(io, "EllipticalFracture{", T, "}(center = ", Tuple(fracture.center),
        ", radii = ", fracture.radii, ", segments = ", fracture.segments,
        ", aperture = ", fracture.properties.aperture)
    show_fracture_mesh_options(io, fracture.mesh_options)
    print(io, ")")
end

function Base.show(io::IO, fracture::PolygonFracture{T}) where T
    print(io, "PolygonFracture{", T, "}(")
    show_fracture_count(io, length(fracture.vertices), "vertex"; plural = "vertices")
    print(io, ", aperture = ", fracture.properties.aperture)
    show_fracture_mesh_options(io, fracture.mesh_options)
    print(io, ")")
end

function show_fracture_mesh_options(io, options)
    if !isnothing(options.resolution)
        print(io, ", resolution = ", options.resolution)
    elseif !isnothing(options.cell_length)
        print(io, ", cell_length = ", options.cell_length)
    end
end
