export reservoir_mesh
# Convenience constructors where you do not have to pass the mesh type
# explicitly. The mesh type is inferred from the argument type.

"""
    mesh = reservoir_mesh(arg...; kwarg...)

This function has two main modes of usage:

 1. To provide a convenient interface for constructing reservoir meshes without
    having to specify the mesh type explicitly. Unless explicitly noted, the
    mesh will always be returned as an (`UnstructuredMesh`)[@ref] from `Jutul`
    with `z_is_depth = true` to match the `JutulDarcy` mesh convention where the
    third coordinate is oriented downwards (primarily for plotting purposes).
 2. To provide a way to extract the reservoir mesh from a case or model that
    contains a reservoir model.

# Useful features of the `UnstructuredMesh` type

- The (`UnstructuredMesh`)[@ref] type uses a node list and several connectivity
  lists to represent the mesh. This means that we can manipulate the field
  `node_points` which will typically be, as an example, a `Vector{SVector{3,
  Float64}}` for 3D meshes with double precision coordinates. If we want to
  change the coordinates of the mesh, we can apply transformations to the
  `node_points` directly, and all cells will remain connected as before provided
  that the transformations are sufficiently smooth to not cause
  self-intersection of cells.
- The function (`extract_submesh`)[@ref] can be used to extract a submesh to
  remove inactive cells. Alternatively, this can be done in
  (`reservoir_domain`)[@ref] by passing an active cell mask or by setting the
  passed `porosity` to zero.
- The function (`grid_dims_ijk`)[@ref] can be used to get the dimensions of the
  mesh in the i, j, k directions when the mesh at some started out as a
  Cartesian mesh (e.g. cartesian or GRDECL variants).
- Individual cells can be looked up by (`cell_ijk`)[@ref] to get the i, j, k
  indices of a cell in a structured grid. This is useful for setting up boundary
  conditions and wells.
"""
function reservoir_mesh

end

function reservoir_mesh(x::Symbol, arg...; kwarg...)
    return reservoir_mesh(Val(x), arg...; kwarg...)
end

function reservoir_mesh(x::AbstractString, arg...; kwarg...)
    return reservoir_mesh(Symbol(x), arg...; kwarg...)
end

"""
    mesh = reservoir_mesh(nx = 3, ny = 5, nz = 10, Lx = 30.0, Ly = 40.0, Lz = 10.0)

Construct a cartesian mesh with `nx` by `ny` by `nz` cells and physical
dimensions `Lx` by `Ly` by `Lz`. The keywords for `nx` and `ny` must be
specified. The keyword for `nz` is optional and defaults to 1, and the keywords
for `Lx`, `Ly` and `Lz` are optional and default to 1.0 meter. This routine will
always provide a 3D mesh.
"""
function reservoir_mesh(::Val{:cartesian}; nx, ny, nz = 1, Lx = 1.0, Ly = 1.0, Lz = 1.0, kwarg...)
    return reservoir_mesh((nx, ny, nz), (Lx, Ly, Lz); kwarg...)
end

function reservoir_mesh(arg...; kwarg...)
    is_grdecl = :coord in keys(kwarg)
    if is_grdecl
        m = reservoir_mesh(Val(:grdecl), arg...; kwarg...)
    else
        m = reservoir_mesh(Val(:cartesian), arg...; kwarg...)
    end
    return m
end

"""
    mesh = reservoir_mesh(dims::Tuple, lengths::Tuple = ones(length(dims)); kwarg...)
    mesh = reservoir_mesh("cartesian", dims, lengths = ones(length(dims)); kwarg...)
    mesh = reservoir_mesh(:cartesian, dims, lengths = ones(length(dims)); kwarg...)

Construct a cartesian mesh with dimensions `dims` and physical dimensions
`lengths`. The keyword `lengths` is optional and defaults to a tuple of ones
with the same length as `dims`. This routine can produce a 2D mesh if only two
entries are passed.

Optionally an `actnum` keyword can be passed to specify active cells. `actnum`
can either be a 3D array with the same dimensions as `dims` or a vector with the
same length as the number of cells in the mesh. In either case, cells with a
value of zero in `actnum` will be removed from the mesh using `extract_submesh`.
"""
function reservoir_mesh(::Val{:cartesian}, dim::Tuple, arg...; origin = nothing, actnum = missing, kwarg...)
    cm = CartesianMesh(dim, arg...; origin = origin)
    m = UnstructuredMesh(cm; z_is_depth = true, kwarg...)
    if !ismissing(actnum)
        if actnum isa AbstractMatrix
            size(actnum) == grid_dims_ijk(m) || error("Size of actnum must match grid dimensions of the mesh.")
            keep = Int[]
            for c in 1:number_of_cells(m)
                i, j, k = cell_ijk(m, c)
                if actnum[i, j, k] != 0
                    push!(keep, c)
                end
            end
        else
            length(actnum) == number_of_cells(m) || error("Length of actnum must match number of cells in the mesh.")
            keep = findall(i -> i > 0, vec(actnum))
        end
        m = extract_submesh(m, keep, z_is_depth = true)
    end
    return m
end

function reservoir_mesh(::Val{:cartesian}, d::Int; kwarg...)
    # Convenience version for 1D mesh.
    return reservoir_mesh(Val(:cartesian), (d,); kwarg...)
end

function reservoir_mesh(::Val{:cartesian}, d::Int, sz; kwarg...)
    # Convenience version for 1D mesh.
    return reservoir_mesh(Val(:cartesian), (d,), (sz, ); kwarg...)
end

"""
    mesh = reservoir_mesh(grid::AbstractDict; kwarg...)

Construct a mesh from a grid section (similar to a GRDECL file) that contains
either COORD/ZCORN or TOPS/DX/DY/DZ.
"""
function reservoir_mesh(d::AbstractDict; kwarg...)
    return reservoir_mesh(Val(:grdecl), d; kwarg...)
end

function reservoir_mesh(::Val{:grdecl}, d::AbstractDict; kwarg...)
    return mesh_from_grid_section(d, kwarg...)
end

# Convenience functions to get the mesh from a case or model
"""
    mesh = reservoir_mesh(case::JutulCase)

Get the reservoir mesh from a case containing a reservoir model.
"""
function reservoir_mesh(case::JutulCase)
    return reservoir_mesh(case.model)
end

"""
    mesh = reservoir_mesh(model::MultiModel)

Get the reservoir mesh from a model containing a reservoir and wells.
"""
function reservoir_mesh(model::MultiModel)
    rmodel = reservoir_model(model)
    return reservoir_mesh(rmodel)
end

"""
    mesh = reservoir_mesh(model::SimulationModel)

Get the reservoir mesh from a reservoir model.
"""
function reservoir_mesh(model::SimulationModel)
    domain = reservoir_domain(model)
    return reservoir_mesh(domain)
end

"""
    mesh = reservoir_mesh(res_domain::DataDomain)

Get the reservoir mesh from a reservoir domain.
"""
function reservoir_mesh(dd::DataDomain)
    return physical_representation(dd)
end