
# Documentation from Jutul.jl {#Documentation-from-Jutul.jl}

`JutulDarcy.jl` builds upon `Jutul.jl`, which takes care of the heavy lifting in terms of meshes, discretizations and solvers. You can use `JutulDarcy.jl` without knowing the inner workings of `Jutul.jl`, but if you want to dive under the hood the [Jutul.jl manual](https://sintefmath.github.io/Jutul.jl/dev/) and [Jutul.jl docstrings](https://sintefmath.github.io/Jutul.jl/dev/docstrings/) may be useful.

We include the docstrings here for your convenience:
<details class='jldocstring custom-block' open>
<summary><a id='Jutul.AMGPreconditioner' href='#Jutul.AMGPreconditioner'><span class="jlbinding">Jutul.AMGPreconditioner</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



AMG on CPU (Julia native)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/linsolve/precond/amg.jl#L2-L4" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.BlockMajorLayout' href='#Jutul.BlockMajorLayout'><span class="jlbinding">Jutul.BlockMajorLayout</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Same as [`EntityMajorLayout`](/ref/jutul#Jutul.EntityMajorLayout), but the system is a sparse matrix where each entry is a small dense matrix.

For a test system with primary variables P, S and equations E1, E2 and two cells this will give a diagonal of length 2: [(∂E1/∂p)₁ (∂E1/∂S)₁ ; (∂E2/∂p)₁ (∂E2/∂S)₁] [(∂E1/∂p)₂ (∂E1/∂S)₂ ; (∂E2/∂p)₂ (∂E2/∂S)₂]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L131-L139" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.BoundaryFaces' href='#Jutul.BoundaryFaces'><span class="jlbinding">Jutul.BoundaryFaces</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Entity for faces on the boundary (faces that are only connected to a single [`Cells`](/ref/jutul#Jutul.Cells))


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L463-L465" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.CartesianMesh' href='#Jutul.CartesianMesh'><span class="jlbinding">Jutul.CartesianMesh</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CartesianMesh(dims, [Δ, [origin]])
```


Create a Cartesian mesh with dimensions specified by the `Tuple` `dims`.

**Arguments**
- `dims::Tuple`: Number of grid cells in each direction. For example, `(nx, ny)` will give a 2D grids with `nx` cells in the x-direction.
  
- `Δ::Tuple=Tuple(ones(length(dims)))`: Equal length to `dims`. First option: A
  

`Tuple` of scalars where each entry is the length of each cell in that direction. For example, specifying `(Δx, Δy) for a uniform grid with each grid cell having area of`Δx*Δy`. Second option: A`Tuple` of vectors where each entry contains the cell sizes in the direction.
- `origin=zeros(length(dims))`: The origin of the first corner in the grid.
  

**Examples**

Generate a uniform 3D mesh that discretizes a domain of 2 by 3 by 5 units with 3 by 5 by 2 cells:

```julia
julia> CartesianMesh((3, 5, 2), (2.0, 3.0, 5.0))
CartesianMesh (3D) with 3x5x2=30 cells
```


Generate a non-uniform 2D mesh:

```julia
julia> CartesianMesh((2, 3), ([1.0, 2.0], [0.1, 3.0, 2.5]))
CartesianMesh (2D) with 2x3x1=6 cells
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/cart.jl#L3-L29" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.Cells' href='#Jutul.Cells'><span class="jlbinding">Jutul.Cells</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Entity for Cells (closed volumes with averaged properties for a finite-volume solver)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L453-L455" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.CoarseMesh-Tuple{Any, Any}' href='#Jutul.CoarseMesh-Tuple{Any, Any}'><span class="jlbinding">Jutul.CoarseMesh</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
CoarseMesh(G::JutulMesh, p)
```


Construct a coarse mesh from a given `JutulMesh` that can be converted to an `UnstructuredMesh` instance. The second argument `p` should be a partition Vector with one entry per cell in the original grid that assigns that cell to a coarse block. Should be one-indexed and the numbering should be sequential and contain at least one fine cell for each coarse index. This is tested by the function.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/coarse.jl#L14-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.CompactAutoDiffCache' href='#Jutul.CompactAutoDiffCache'><span class="jlbinding">Jutul.CompactAutoDiffCache</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Cache that holds an AD vector/matrix together with their positions.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L663-L665" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.DataDomain-Tuple{JutulDomain}' href='#Jutul.DataDomain-Tuple{JutulDomain}'><span class="jlbinding">Jutul.DataDomain</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
DataDomain(domain::JutulDomain; property1 = p1, property2 = p2, ...)
```


A wrapper around a domain that allows for storing of entity-associated data.

Example:

```julia
# Grid with 6 cells and 7 interior faces
g = CartesianMesh((2, 3))
d = DataDomain(g)
d[:cell_vec] = rand(6) #ok, same as:
d[:cell_vec, Cells()] = rand(6) #ok
d[:cell_vec, Faces()] = rand(6) #not ok!
d[:face_vec, Faces()] = rand(7) #ok!
# Can also add general arrays if last dimension == entity dimension
d[:cell_vec, Cells()] = rand(10, 3, 6) #ok
# Can add general data too, but needs to be specified
d[:not_on_face_or_cell, nothing] = rand(3) # also ok
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/domains.jl#L119-L138" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.DefaultContext' href='#Jutul.DefaultContext'><span class="jlbinding">Jutul.DefaultContext</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Default context


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/contexts/default.jl#L2" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.DiscretizedDomain' href='#Jutul.DiscretizedDomain'><span class="jlbinding">Jutul.DiscretizedDomain</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DiscretizedDomain(domain, disc = nothing)
```


A type for a discretized domain of some other domain or mesh. May contain one or more discretizations as-needed to write equations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/domains.jl#L40-L45" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.EntityMajorLayout' href='#Jutul.EntityMajorLayout'><span class="jlbinding">Jutul.EntityMajorLayout</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Equations are grouped by entity, listing all equations and derivatives for entity 1 before proceeding to entity 2 etc.

For a test system with primary variables P, S and equations E1, E2 and two cells this will give the following ordering on the diagonal: (∂E1/∂p)₁, (∂E2/∂S)₁, (∂E1/∂p)₂, (∂E2/∂S)₂


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L115-L122" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.EquationMajorLayout' href='#Jutul.EquationMajorLayout'><span class="jlbinding">Jutul.EquationMajorLayout</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Equations are stored sequentially in rows, derivatives of same type in columns:

For a test system with primary variables P, S and equations E1, E2 and two cells this will give the following ordering on the diagonal: (∂E1/∂p)₁, (∂E1/∂p)₂, (∂E2/∂S)₁, (∂E2/∂S)₂


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L102-L108" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.EquationSet' href='#Jutul.EquationSet'><span class="jlbinding">Jutul.EquationSet</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Set of a variable where equations are defined


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/dd/dd.jl#L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.FaceMap' href='#Jutul.FaceMap'><span class="jlbinding">Jutul.FaceMap</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Struct that contains mappings for a set of faces that are made up of nodes and are part of cells.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/unstructured/types.jl#L1-L4" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.Faces' href='#Jutul.Faces'><span class="jlbinding">Jutul.Faces</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Entity for Faces (intersection between pairs of [`Cells`](/ref/jutul#Jutul.Cells))


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L458-L460" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.FlowDiscretization' href='#Jutul.FlowDiscretization'><span class="jlbinding">Jutul.FlowDiscretization</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Discretization of kgradp + upwind


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L762" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.FractionVariables' href='#Jutul.FractionVariables'><span class="jlbinding">Jutul.FractionVariables</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for fraction variables (vector variables that sum up to unity over each entity).

By default, these are limited to the [0, 1] range through [`maximum_value`](/ref/jutul#Jutul.maximum_value-Tuple{JutulVariables}) and [`minimum_value`](/ref/jutul#Jutul.minimum_value-Tuple{JutulVariables}) default implementations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L69-L75" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.GenericKrylov' href='#Jutul.GenericKrylov'><span class="jlbinding">Jutul.GenericKrylov</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GenericKrylov(solver = :gmres; preconditioner = nothing; <kwarg>)
```


Solver that wraps `Krylov.jl` with support for preconditioning.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/linsolve/krylov.jl#L29-L33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.GlobalSet' href='#Jutul.GlobalSet'><span class="jlbinding">Jutul.GlobalSet</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



The global set of variables


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/dd/dd.jl#L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.GroupWisePreconditioner' href='#Jutul.GroupWisePreconditioner'><span class="jlbinding">Jutul.GroupWisePreconditioner</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Multi-model preconditioners


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/linsolve/precond/various.jl#L58-L60" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.HalfFaces' href='#Jutul.HalfFaces'><span class="jlbinding">Jutul.HalfFaces</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Entity for half-faces (face associated with a single [`Cells`](/ref/jutul#Jutul.Cells))


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L468-L470" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.HelperSimulator-Union{Tuple{E}, Tuple{M}, Tuple{M, Any}} where {M, E}' href='#Jutul.HelperSimulator-Union{Tuple{E}, Tuple{M}, Tuple{M, Any}} where {M, E}'><span class="jlbinding">Jutul.HelperSimulator</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
HelperSimulator(model::M, T = Float64; state0 = setup_state(model), executor::E = Jutul.default_executor()) where {M, E}
```


Construct a helper simulator that can be used to compute the residuals and/or accumulation terms for a given type T. Useful for coupling Jutul to other solvers and types of automatic differentiation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/helper.jl#L9-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ILUZeroPreconditioner' href='#Jutul.ILUZeroPreconditioner'><span class="jlbinding">Jutul.ILUZeroPreconditioner</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



ILU(0) preconditioner on CPU


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/linsolve/precond/ilu.jl#L1-L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.IndirectionMap' href='#Jutul.IndirectionMap'><span class="jlbinding">Jutul.IndirectionMap</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



IndirectionMap(vals::Vector{V}, pos::Vector{Int}) where V

Create a indirection map that encodes a variable length dense vector.

`pos` is assumed to be a Vector{Int} of length n+1 where n is the number of dense vectors that is encoded. The `vals` array holds the entries for vector i in the range `pos[i]:(pos[i+1]-1)` for fast lookup. Indexing into the indirection map with index `k` will give a view into the values for vector `k`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L1073-L1082" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.IndirectionMap-Union{Tuple{Array{Vector{T}, 1}}, Tuple{T}} where T' href='#Jutul.IndirectionMap-Union{Tuple{Array{Vector{T}, 1}}, Tuple{T}} where T'><span class="jlbinding">Jutul.IndirectionMap</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
imap = IndirectionMap(vec_of_vec::Vector{V}) where V
```


Create indirection map for a variable length dense vector that is represented as a Vector of Vectors.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L1094-L1099" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JacobiPreconditioner' href='#Jutul.JacobiPreconditioner'><span class="jlbinding">Jutul.JacobiPreconditioner</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Damped Jacobi preconditioner on CPU


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/linsolve/precond/jacobi.jl#L2-L4" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulAutoDiffCache' href='#Jutul.JutulAutoDiffCache'><span class="jlbinding">Jutul.JutulAutoDiffCache</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



An AutoDiffCache is a type that holds both a set of AD values and a map into some global Jacobian.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L658-L661" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulCase' href='#Jutul.JutulCase'><span class="jlbinding">Jutul.JutulCase</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
JutulCase(model, dt = [1.0], forces = setup_forces(model); state0 = nothing, parameters = nothing, kwarg...)
```


Set up a structure that holds the complete specification of a simulation case.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L840-L844" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulConfig' href='#Jutul.JutulConfig'><span class="jlbinding">Jutul.JutulConfig</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



JutulConfig(name = nothing)

A configuration object that acts like a `Dict{Symbol,Any}` but contains additional data to limit the valid keys and values to those added by [`add_option!`](/ref/jutul#Jutul.add_option!)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/types.jl#L117-L122" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulContext' href='#Jutul.JutulContext'><span class="jlbinding">Jutul.JutulContext</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for the context Jutul should execute in (matrix formats, memory allocation, etc.)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L83-L85" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulDiscretization' href='#Jutul.JutulDiscretization'><span class="jlbinding">Jutul.JutulDiscretization</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
d = disc(i, Cells())
```


Ask discretization for entry `i` when discretizing some equation on the  chosen entity (e.g. [`Cells`](/ref/jutul#Jutul.Cells))


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L26-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulDiscretization-2' href='#Jutul.JutulDiscretization-2'><span class="jlbinding">Jutul.JutulDiscretization</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for a Jutul discretization


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L21-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulDomain' href='#Jutul.JutulDomain'><span class="jlbinding">Jutul.JutulDomain</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for domains where equations can be defined


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/domains.jl#L1-L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulEntity' href='#Jutul.JutulEntity'><span class="jlbinding">Jutul.JutulEntity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Super-type for all entities where [`JutulVariables`](/ref/jutul#Jutul.JutulVariables) can be defined.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L448-L450" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulEquation' href='#Jutul.JutulEquation'><span class="jlbinding">Jutul.JutulEquation</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for all residual equations


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L229-L231" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulForce' href='#Jutul.JutulForce'><span class="jlbinding">Jutul.JutulForce</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for driving forces


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L78-L80" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulMatrixLayout' href='#Jutul.JutulMatrixLayout'><span class="jlbinding">Jutul.JutulMatrixLayout</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for matrix layouts. A layout determines how primary variables and equations are ordered in a sparse matrix representation. Note that this is different from the matrix format itself as it concerns the ordering itself: For example, if all equations for a single cell come in sequence, or if a single equation is given for all entities before the next equation is written.

Different layouts does not change the solution of the system, but different linear solvers support different layouts.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L91-L100" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulMesh' href='#Jutul.JutulMesh'><span class="jlbinding">Jutul.JutulMesh</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



A mesh is a type of domain that has been discretized. Abstract subtype.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L442-L444" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulSystem' href='#Jutul.JutulSystem'><span class="jlbinding">Jutul.JutulSystem</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for the physical system to be solved.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L16-L18" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.JutulVariables' href='#Jutul.JutulVariables'><span class="jlbinding">Jutul.JutulVariables</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for all variables in Jutul.

A variable is associated with a [`JutulEntity`](/ref/jutul#Jutul.JutulEntity) through the [`associated_entity`](/ref/jutul#Jutul.associated_entity-Tuple{JutulEquation}) function. A variable is local to that entity, and cannot depend on other entities. Variables are used by models to define:
- primary variables: Sometimes referred to as degrees of freedom, primary unknowns or solution variables
  
- parameters: Static quantities that impact the solution
  
- secondary variables: Can be computed from a combination of other primary and secondary variables and parameters.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L45-L56" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.LUPreconditioner' href='#Jutul.LUPreconditioner'><span class="jlbinding">Jutul.LUPreconditioner</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Full LU factorization as preconditioner (intended for smaller subsystems)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/linsolve/precond/various.jl#L8-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.LUSolver' href='#Jutul.LUSolver'><span class="jlbinding">Jutul.LUSolver</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
LUSolver(; reuse_memory = true, check = true, max_size = 50000)
```


Direct solver that calls `lu` directly. Direct solvers are highly accurate, but are costly in terms of memory usage and execution speed for larger systems.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/linsolve/scalar_cpu.jl#L2-L7" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.LimitByFailedTimestepSelector' href='#Jutul.LimitByFailedTimestepSelector'><span class="jlbinding">Jutul.LimitByFailedTimestepSelector</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
sel = LimitByFailedTimestepSelector(num = 10, factor = 0.9)
```


Limit the timestep by the shortest of `num` failed timesteps, reducing the timestep by `factor` multiplied by the shortest failed timestep. If no time-steps failed during the last `num` steps, the timestep is not changed.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/timesteps.jl#L151-L157" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.MRSTWrapMesh' href='#Jutul.MRSTWrapMesh'><span class="jlbinding">Jutul.MRSTWrapMesh</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MRSTWrapMesh(G, N = nothing)
```


Mesh that adapts an exported MRST mesh to the Jutul interface. `G` is assumed to be read directly from file using `MAT.matread`. The raw exported grid can be found under the `data` field.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/mrst.jl#L9-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.MultiModel' href='#Jutul.MultiModel'><span class="jlbinding">Jutul.MultiModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MultiModel(models)
MultiModel(models, :SomeLabel)
```


A model variant that is made up of many named submodels, each a fully realized [`SimulationModel`](/ref/jutul#Jutul.SimulationModel-Tuple{Any,%20Any}).

`models` should be a `NamedTuple` or `Dict{Symbol, JutulModel}`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L954-L961" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.NoEntity' href='#Jutul.NoEntity'><span class="jlbinding">Jutul.NoEntity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



An entity for something that isn&#39;t associated with an entity


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L478-L480" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.Nodes' href='#Jutul.Nodes'><span class="jlbinding">Jutul.Nodes</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Entity for Nodes (intersection between multiple [`Faces`](/ref/jutul#Jutul.Faces))


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L473-L475" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ParallelCSRContext' href='#Jutul.ParallelCSRContext'><span class="jlbinding">Jutul.ParallelCSRContext</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



A context that uses a CSR sparse matrix format together with threads. Experimental.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/contexts/csr.jl#L2" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.SPAI0Preconditioner' href='#Jutul.SPAI0Preconditioner'><span class="jlbinding">Jutul.SPAI0Preconditioner</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Sparse Approximate Inverse preconditioner of lowest order – SPAI(0)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/linsolve/precond/spai.jl#L1-L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.SPU' href='#Jutul.SPU'><span class="jlbinding">Jutul.SPU</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Single-point upwinding.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/conservation/flux.jl#L31-L33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ScalarVariable' href='#Jutul.ScalarVariable'><span class="jlbinding">Jutul.ScalarVariable</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for scalar variables (one entry per entity, e.g. pressure or temperature in each cell of a model)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L59-L62" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.SimulationModel-Tuple{Any, Any}' href='#Jutul.SimulationModel-Tuple{Any, Any}'><span class="jlbinding">Jutul.SimulationModel</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
SimulationModel(domain, system; <kwarg>)
```


Instantiate a model for a given `system` discretized on the `domain`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L259-L263" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.SimulationModel-Tuple{JutulMesh, Any}' href='#Jutul.SimulationModel-Tuple{JutulMesh, Any}'><span class="jlbinding">Jutul.SimulationModel</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
SimulationModel(g::JutulMesh, system; discretization = nothing, kwarg...)
```


Type that defines a simulation model - everything needed to solve the discrete equations.

The minimal setup requires a [`JutulMesh`](/ref/jutul#Jutul.JutulMesh) that defines topology together with a [`JutulSystem`](/ref/jutul#Jutul.JutulSystem) that imposes physical laws.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L486-L494" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.Simulator-Tuple{Any}' href='#Jutul.Simulator-Tuple{Any}'><span class="jlbinding">Jutul.Simulator</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
Simulator(model; <kwarg>)
```


Set up a simulator object for a `model` that can be used by [`simulate!`](/ref/jutul#Jutul.simulate!-Tuple{JutulSimulator,%20AbstractVector}). To avoid manually instantiating the simulator, the non-mutating [`simulate`](/ref/jutul#Jutul.simulate-Tuple{Any,%20JutulModel,%20AbstractVector}) interface can be used instead.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/types.jl#L19-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.SparsityTracingWrapper-Union{Tuple{AbstractArray{T, N}}, Tuple{N}, Tuple{T}} where {T, N}' href='#Jutul.SparsityTracingWrapper-Union{Tuple{AbstractArray{T, N}}, Tuple{N}, Tuple{T}} where {T, N}'><span class="jlbinding">Jutul.SparsityTracingWrapper</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
SparsityTracingWrapper(x::AbstractArray{T, N}) where {T, N}
```


Create a sparsity tracing wrapper for a numeric array. This wrapped array produces outputs that have the same value as the wrapped type, but contains a SparsityTracing seeded value with seed equal to the column index (if matrix) or linear index (if vector).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/sparsity.jl#L22-L29" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.TPFA' href='#Jutul.TPFA'><span class="jlbinding">Jutul.TPFA</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Two-point flux approximation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/conservation/flux.jl#L12-L14" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.TwoPointFiniteVolumeGeometry' href='#Jutul.TwoPointFiniteVolumeGeometry'><span class="jlbinding">Jutul.TwoPointFiniteVolumeGeometry</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TwoPointFiniteVolumeGeometry(neighbors, areas, volumes, normals, cell_centers, face_centers)
```


Store two-point geometry information for a given list of `neighbors` specified as a `2` by `n` matrix where `n` is the number of faces such that face `i` connectes cells `N[1, i]` and `N[2, i]`.

The two-point finite-volume geometry contains the minimal set of geometry information required to compute standard finite-volume discretizations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/meshes.jl#L19-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.UnstructuredMesh-Tuple{CartesianMesh}' href='#Jutul.UnstructuredMesh-Tuple{CartesianMesh}'><span class="jlbinding">Jutul.UnstructuredMesh</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
UnstructuredMesh(g::CartesianMesh)
```


Convert `CartesianMesh` instance to unstructured grid. Note that the mesh must be 2D and 3D for a 1-to-1 conversion. 1D meshes are implicitly converted to 2D.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/unstructured/types.jl#L289-L294" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.VariableSet' href='#Jutul.VariableSet'><span class="jlbinding">Jutul.VariableSet</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Set of a variable where variables are defined


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/dd/dd.jl#L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.VectorVariables' href='#Jutul.VectorVariables'><span class="jlbinding">Jutul.VectorVariables</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for vector variables (more than one entry per entity, for example saturations or displacements)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L64-L67" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.absolute_increment_limit-Tuple{JutulVariables}' href='#Jutul.absolute_increment_limit-Tuple{JutulVariables}'><span class="jlbinding">Jutul.absolute_increment_limit</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Absolute allowable change for variable during a nonlinear update.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L86-L88" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.add_option!' href='#Jutul.add_option!'><span class="jlbinding">Jutul.add_option!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
add_option!(opts::JutulConfig, :my_cool_option, 3, "My option has this brief description")
```


Add an option to existing [`JutulConfig`](/ref/jutul#Jutul.JutulConfig) structure. Additional currently undocumented keyword arguments can be used to restrict valid types and values.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/config.jl#L3-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.align_to_jacobian!-NTuple{4, Any}' href='#Jutul.align_to_jacobian!-NTuple{4, Any}'><span class="jlbinding">Jutul.align_to_jacobian!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Update an equation so that it knows where to store its derivatives in the Jacobian representation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L422-L425" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.align_to_jacobian!-Tuple{ConservationLawTPFAStorage, ConservationLaw, Any, Any, Cells}' href='#Jutul.align_to_jacobian!-Tuple{ConservationLawTPFAStorage, ConservationLaw, Any, Any, Cells}'><span class="jlbinding">Jutul.align_to_jacobian!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Update positions of law&#39;s derivatives in global Jacobian


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/conservation/conservation.jl#L82" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.allocate_array_ad-Tuple{AbstractMatrix}' href='#Jutul.allocate_array_ad-Tuple{AbstractMatrix}'><span class="jlbinding">Jutul.allocate_array_ad</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
allocate_array_ad(v::AbstractMatrix, ...)
```


Convert matrix to AD matrix.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L346-L349" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.allocate_array_ad-Tuple{AbstractVector}' href='#Jutul.allocate_array_ad-Tuple{AbstractVector}'><span class="jlbinding">Jutul.allocate_array_ad</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
allocate_array_ad(v::AbstractVector, ...)
```


Convert vector to AD vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L336-L339" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.allocate_array_ad-Union{Tuple{Vararg{R}}, Tuple{R}} where R<:Integer' href='#Jutul.allocate_array_ad-Union{Tuple{Vararg{R}}, Tuple{R}} where R<:Integer'><span class="jlbinding">Jutul.allocate_array_ad</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
allocate_array_ad(n[, m]; <keyword arguments>)
```


Allocate vector or matrix as AD with optionally provided context and a specified non-zero on the diagonal.

**Arguments**
- `n::Integer`: number of entries in vector, or number of rows if `m` is given.
  
- `m::Integer`: number of rows (optional)
  

**Keyword arguments**
- `npartials = 1`: Number of partials derivatives to allocate for each element
  
- `diag_pos = nothing`: Indices of where to put entities on the diagonal (if any)
  

Other keyword arguments are passed onto `get_ad_entity_scalar`.

**Examples:**

Allocate a vector with a single partial:

```julia
julia> allocate_array_ad(2)
2-element Vector{ForwardDiff.Dual{nothing, Float64, 1}}:
 Dual{nothing}(0.0,0.0)
 Dual{nothing}(0.0,0.0)
```


Allocate a vector with two partials, and set the first to one:

```julia
julia> allocate_array_ad(2, diag_pos = 1, npartials = 2)
2-element Vector{ForwardDiff.Dual{nothing, Float64, 2}}:
 Dual{nothing}(0.0,1.0,0.0)
 Dual{nothing}(0.0,1.0,0.0)
```


Set up a matrix with two partials, where the first column has partials [1, 0] and the second [0, 1]:

```julia
julia> allocate_array_ad(2, 2, diag_pos = [1, 2], npartials = 2)
2×2 Matrix{ForwardDiff.Dual{nothing, Float64, 2}}:
 Dual{nothing}(0.0,1.0,0.0)  Dual{nothing}(0.0,1.0,0.0)
 Dual{nothing}(0.0,0.0,1.0)  Dual{nothing}(0.0,0.0,1.0)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L278-L316" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.apply!-Tuple{Any, TrivialPreconditioner, Any, Vararg{Any}}' href='#Jutul.apply!-Tuple{Any, TrivialPreconditioner, Any, Vararg{Any}}'><span class="jlbinding">Jutul.apply!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Trivial / identity preconditioner with size for use in subsystems.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/linsolve/precond/various.jl#L37-L39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.apply_forces!-NTuple{4, Any}' href='#Jutul.apply_forces!-NTuple{4, Any}'><span class="jlbinding">Jutul.apply_forces!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Apply a set of forces to all equations. Equations that don&#39;t support a given force will just ignore them, thanks to the power of multiple dispatch.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L853-L856" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.apply_forces_to_equation!-NTuple{7, Any}' href='#Jutul.apply_forces_to_equation!-NTuple{7, Any}'><span class="jlbinding">Jutul.apply_forces_to_equation!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
apply_forces_to_equation!(diag_part, storage, model, eq, eq_s, force, time)
```


Update an equation with the effect of a force. The default behavior for any force we do not know about is to assume that the force does not impact this particular equation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L570-L576" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.as_value-Tuple{AbstractArray}' href='#Jutul.as_value-Tuple{AbstractArray}'><span class="jlbinding">Jutul.as_value</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Create a mapped array that produces only the values when indexed.

Only useful for AD arrays, otherwise it does nothing.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L437-L441" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.associated_entity-Tuple{JutulEquation}' href='#Jutul.associated_entity-Tuple{JutulEquation}'><span class="jlbinding">Jutul.associated_entity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Return the domain entity the equation is associated with


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L252-L254" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.associated_entity-Tuple{JutulVariables}' href='#Jutul.associated_entity-Tuple{JutulVariables}'><span class="jlbinding">Jutul.associated_entity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



The entity a variable is associated with, and can hold partial derivatives with respect to.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L9-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.cartesian_partition' href='#Jutul.cartesian_partition'><span class="jlbinding">Jutul.cartesian_partition</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
cartesian_partition(D::DataDomain, coarse_dims, ptype = :centroids)
```


Perform Cartesian partition in centroid or IJK space.

**Examples**

Generate 5x5x5 partition based on points (suitable for unstructured meshes)

```julia
mesh = CartesianMesh((50, 50, 1))
p_pts = Jutul.cartesian_partition(domain, (5, 5, 5), :centroids)
```


Generate 5x5x5 partition based on IJK indices (suitable for logically Cartesian meshes)

```julia
mesh = CartesianMesh((50, 50, 1))
p_ijk = Jutul.cartesian_partition(domain, (5, 5, 5), :ijk)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/partitioning.jl#L163-L184" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.cell_dims-Tuple{Any, Any}' href='#Jutul.cell_dims-Tuple{Any, Any}'><span class="jlbinding">Jutul.cell_dims</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cell_dims(g, pos)::Tuple
```


Get physical cell dimensions of cell with index `pos` for grid `g`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/cart.jl#L130-L134" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.cell_index-Tuple{Any, Tuple}' href='#Jutul.cell_index-Tuple{Any, Tuple}'><span class="jlbinding">Jutul.cell_index</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cell_index(g, pos)
```


Get linear (scalar) index of mesh cell from provided IJK tuple `pos`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/cart.jl#L108-L112" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.cells_inside_bounding_box-Tuple{UnstructuredMesh, Any, Any}' href='#Jutul.cells_inside_bounding_box-Tuple{UnstructuredMesh, Any, Any}'><span class="jlbinding">Jutul.cells_inside_bounding_box</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cells_inside_bounding_box(G::UnstructuredMesh, low_bb, high_bb; algorithm = :box, atol = 0.01)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/trajectories.jl#L207-L211" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.check_amgcl_availability-Tuple{}' href='#Jutul.check_amgcl_availability-Tuple{}'><span class="jlbinding">Jutul.check_amgcl_availability</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
check_amgcl_availability(; throw = true)
```


Check if AMGCLWrap extension is available. If `throw=true` this wil be an error, otherwise a Boolean indicating if the extension is available will be returned.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ext/amgclwrap_ext.jl#L20-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.check_plotting_availability-Tuple{}' href='#Jutul.check_plotting_availability-Tuple{}'><span class="jlbinding">Jutul.check_plotting_availability</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
check_plotting_availability(; throw = true, interactive = false)
```


Check if plotting through at least one `Makie` backend is available in the Julia session (after package has been loaded by for example `using GLMakie`). The argument `throw` can be used to control if this function acts as a programmatic check (`throw=false`) there the return value indicates availability, or if an error message is to be printed telling the user how to get plotting working (`throw=true`)

An additional check for specifically `interactive` plots can also be added.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ext/makie_ext.jl#L141-L152" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.compress_timesteps' href='#Jutul.compress_timesteps'><span class="jlbinding">Jutul.compress_timesteps</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compress_timesteps(timesteps, forces = nothing; max_step = Inf)
```


Compress a set of timesteps and forces to the largest possible steps that still covers the same interval and changes forces at exactly the same points in time, while being limited to a maximum size of `max_step`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/timesteps.jl#L248-L254" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.compress_timesteps-Tuple{JutulCase}' href='#Jutul.compress_timesteps-Tuple{JutulCase}'><span class="jlbinding">Jutul.compress_timesteps</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
compress_timesteps(case::JutulCase; max_step = Inf)
```


Compress time steps for a Jutul case. See [`compress_timesteps`](/ref/jutul#Jutul.compress_timesteps) for the general case.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/timesteps.jl#L305-L310" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.compute_boundary_trans-Tuple{DataDomain, Any}' href='#Jutul.compute_boundary_trans-Tuple{DataDomain, Any}'><span class="jlbinding">Jutul.compute_boundary_trans</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
compute_boundary_trans(d::DataDomain, perm)
```


Compute the boundary half face transmissibilities for perm. The input `perm` can either be the symbol of some data defined on `Cells()`, a vector of numbers for each cell or a matrix with number of columns equal to the number of cells.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/discretization/finite-volume.jl#L230-L237" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.compute_face_trans-Tuple{DataDomain, Vararg{Any}}' href='#Jutul.compute_face_trans-Tuple{DataDomain, Vararg{Any}}'><span class="jlbinding">Jutul.compute_face_trans</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



compute_face_trans(g::DataDomain, perm)

Compute face trans for the interior faces. The input `perm` can either be the symbol of some data defined on `Cells()`, a vector of numbers for each cell or a matrix with number of columns equal to the number of cells.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/discretization/finite-volume.jl#L213-L219" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.compute_half_face_trans-Tuple{DataDomain, Any}' href='#Jutul.compute_half_face_trans-Tuple{DataDomain, Any}'><span class="jlbinding">Jutul.compute_half_face_trans</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
compute_half_face_trans(g::DataDomain, perm)
```


Compute half-face trans for the interior faces. The input `perm` can either be the symbol of some data defined on `Cells()`, a vector of numbers for each cell or a matrix with number of columns equal to the number of cells.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/discretization/finite-volume.jl#L8-L14" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.convergence_criterion-Tuple{Any, Any, JutulEquation, Any, Any}' href='#Jutul.convergence_criterion-Tuple{Any, Any, JutulEquation, Any, Any}'><span class="jlbinding">Jutul.convergence_criterion</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convergence_criterion(model, storage, eq, eq_s, r; dt = 1)
```


Get the convergence criterion values for a given equation. Can be checked against the corresponding tolerances.

**Arguments**
- `model`: model that generated the current equation.
  
- `storage`: global simulator storage.
  
- `eq::JutulEquation`: equation implementation currently being checked
  
- `eq_s`: storage for `eq` where values are contained.
  
- `r`: the local residual part corresponding to this model, as a matrix with column index equaling entity index
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L581-L592" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.convert_from_si-Tuple{Any, String}' href='#Jutul.convert_from_si-Tuple{Any, String}'><span class="jlbinding">Jutul.convert_from_si</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_from_si(value, unit_name::Union{Symbol, String})
```


Convert `value` from SI representation to the unit in `unit_symbol`.

**Examples**

```julia
julia> convert_from_si(3600.0, :hour) # Get 3600 s represented as hours
1.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/units/interface.jl#L29-L39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.convert_state_ad' href='#Jutul.convert_state_ad'><span class="jlbinding">Jutul.convert_state_ad</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



Convert a state containing variables as arrays of doubles to a state where those arrays contain the same value as Dual types. The dual type is currently taken from ForwardDiff.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L217-L221" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.convert_to_si-Tuple{Any, String}' href='#Jutul.convert_to_si-Tuple{Any, String}'><span class="jlbinding">Jutul.convert_to_si</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_to_si(value, unit_name::String)
```


Convert `value` to SI representation from value in the unit given by `unit_symbol`.

**Examples**

```julia
julia> convert_to_si(1.0, :hour) # Get 1 hour represented as seconds
3600.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/units/interface.jl#L1-L12" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.coord_offset-Tuple{Any, AbstractFloat}' href='#Jutul.coord_offset-Tuple{Any, AbstractFloat}'><span class="jlbinding">Jutul.coord_offset</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Lower corner for one dimension, without any transforms applied


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/cart.jl#L102-L104" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.data_domain_to_parameters_gradient-Tuple{Any, Any}' href='#Jutul.data_domain_to_parameters_gradient-Tuple{Any, Any}'><span class="jlbinding">Jutul.data_domain_to_parameters_gradient</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
data_domain_to_parameters_gradient(model, parameter_gradient; dp_dd = missing, config = nothing)
```


Make a data_domain copy that contains the gradient of some objective with respect to the fields in the data_domain, assuming that the parameters were initialized directly from the data_domain via (`setup_parameters`)[@ref].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/vectorization.jl#L229-L235" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.declare_pattern-Tuple{Any, Any, Any, Any, Vararg{Any}}' href='#Jutul.declare_pattern-Tuple{Any, Any, Any, Any, Vararg{Any}}'><span class="jlbinding">Jutul.declare_pattern</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Give out source, target arrays of equal length for a given equation attached to the given model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L388-L391" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.declare_sparsity' href='#Jutul.declare_sparsity'><span class="jlbinding">Jutul.declare_sparsity</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



Give out I, J arrays of equal length for a given equation attached to the given model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L301-L304" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.degrees_of_freedom_per_entity-Tuple{Any, ScalarVariable}' href='#Jutul.degrees_of_freedom_per_entity-Tuple{Any, ScalarVariable}'><span class="jlbinding">Jutul.degrees_of_freedom_per_entity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Number of independent primary variables / degrees of freedom per computational entity.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L73-L75" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.descalarize_primary_variable!-Tuple{Any, Any, Any, ScalarVariable, Any}' href='#Jutul.descalarize_primary_variable!-Tuple{Any, Any, Any, ScalarVariable, Any}'><span class="jlbinding">Jutul.descalarize_primary_variable!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
descalarize_primary_variable!(dest_array, model, V, var::Jutul.ScalarVariable, index)
```


Descalarize a primary variable, overwriting dest_array at entity `index`. The AD status of entries in `dest_array` will be retained.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/scalarization.jl#L49-L54" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.descalarize_primary_variables!' href='#Jutul.descalarize_primary_variables!'><span class="jlbinding">Jutul.descalarize_primary_variables!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
descalarize_primary_variables!(state, model, V, pvars::NamedTuple = (; pairs(model.primary_variables)...), ind = eachindex(V))
```


Replace valeus in `state` by the scalarized values found in V.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/scalarization.jl#L130-L134" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.dim-Tuple{JutulMesh}' href='#Jutul.dim-Tuple{JutulMesh}'><span class="jlbinding">Jutul.dim</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dim(g)::Integer
```


Get the dimension of a mesh.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/meshes.jl#L136-L140" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.expand_to_ministeps-Tuple{Any, Any}' href='#Jutul.expand_to_ministeps-Tuple{Any, Any}'><span class="jlbinding">Jutul.expand_to_ministeps</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
substates, dt, report_index = expand_to_ministeps(states, reports)
```


Get states and timesteps at the finest stored resolution. Output lengths depend on if `output_substates` option to simulator was enabled.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/utils.jl#L748-L753" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.expand_to_ministeps-Tuple{Jutul.SimResult}' href='#Jutul.expand_to_ministeps-Tuple{Jutul.SimResult}'><span class="jlbinding">Jutul.expand_to_ministeps</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
substates, dt, report_index = expand_to_ministeps(result::SimResult)
```


Get states and timesteps at the finest stored resolution. Output lengths depend on if `output_substates` option to simulator was enabled.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/utils.jl#L737-L742" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.extra_debug_output!-NTuple{6, Any}' href='#Jutul.extra_debug_output!-NTuple{6, Any}'><span class="jlbinding">Jutul.extra_debug_output!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
extra_debug_output!(report, storage, model, config, iteration, dt)
```


Add extra debug output to report during a nonlinear iteration.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/simulator.jl#L756-L760" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.extract_submesh-Tuple{UnstructuredMesh, Any}' href='#Jutul.extract_submesh-Tuple{UnstructuredMesh, Any}'><span class="jlbinding">Jutul.extract_submesh</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
extract_submesh(g::UnstructuredMesh, cells)
```


Extract a subgrid for a given mesh and a iterable of `cells` to keep.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/unstructured/utils.jl#L27-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.extrude_mesh-Tuple{UnstructuredMesh, Int64}' href='#Jutul.extrude_mesh-Tuple{UnstructuredMesh, Int64}'><span class="jlbinding">Jutul.extrude_mesh</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
extrude_mesh(m2d::UnstructuredMesh, nlayers)
extrude_mesh(m2d::UnstructuredMesh, [1, 2, 5, 10])
```


Extrude a 2D mesh into a 3D mesh by adding layers of cells in the z-direction. The number of layers can be specified as an integer or as an array of depths. The depths must be in increasing order.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/extruded.jl#L2-L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.find_enclosing_cell-Union{Tuple{T}, Tuple{D}, Tuple{UnstructuredMesh{D}, StaticArraysCore.SVector{D, T}, Vararg{AbstractArray{StaticArraysCore.SVector{D, T}, 1}, 4}}, Tuple{UnstructuredMesh{D}, StaticArraysCore.SVector{D, T}, AbstractArray{StaticArraysCore.SVector{D, T}, 1}, AbstractArray{StaticArraysCore.SVector{D, T}, 1}, AbstractArray{StaticArraysCore.SVector{D, T}, 1}, AbstractArray{StaticArraysCore.SVector{D, T}, 1}, Any}} where {D, T}' href='#Jutul.find_enclosing_cell-Union{Tuple{T}, Tuple{D}, Tuple{UnstructuredMesh{D}, StaticArraysCore.SVector{D, T}, Vararg{AbstractArray{StaticArraysCore.SVector{D, T}, 1}, 4}}, Tuple{UnstructuredMesh{D}, StaticArraysCore.SVector{D, T}, AbstractArray{StaticArraysCore.SVector{D, T}, 1}, AbstractArray{StaticArraysCore.SVector{D, T}, 1}, AbstractArray{StaticArraysCore.SVector{D, T}, 1}, AbstractArray{StaticArraysCore.SVector{D, T}, 1}, Any}} where {D, T}'><span class="jlbinding">Jutul.find_enclosing_cell</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
find_enclosing_cell(G::UnstructuredMesh{D}, pt::SVector{D, T},
    normals::AbstractVector{SVector{D, T}},
    face_centroids::AbstractVector{SVector{D, T}},
    boundary_normals::AbstractVector{SVector{D, T}},
    boundary_centroids::AbstractVector{SVector{D, T}},
    cells = 1:number_of_cells(G)
) where {D, T}
```


Find enclosing cell of a point. This can be a bit expensive for larger meshes. Recommended to use the more high level `find_enclosing_cells` instead.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/trajectories.jl#L132-L143" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.find_enclosing_cells-Tuple{Any, Any}' href='#Jutul.find_enclosing_cells-Tuple{Any, Any}'><span class="jlbinding">Jutul.find_enclosing_cells</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
find_enclosing_cells(G, traj; geometry = tpfv_geometry(G), n = 25)
```


Find the cell indices of cells in the mesh `G` that are intersected by a given trajectory `traj`. `traj` can be either a matrix with equal number of columns as dimensions in G (i.e. three columns for 3D) or a `Vector` of `SVector` instances with the same length.

The optional argument `geometry` is used to define the centroids and normals used in the calculations. You can precompute this if you need to perform many searches. The keyword argument `n` can be used to set the number of discretizations in each segment.

`use_boundary` is by default set to `false`. If set to true, the boundary faces of cells are treated more rigorously when picking exactly what cells are cut by a trajectory, but this requires that the boundary normals are oriented outwards, which is currently not the case for all meshes from downstream packages.

`limit_box` speeds up the search by limiting the search to the minimal bounding box that contains both the trajectory and the mesh. This can be turned off by passing `false`. There should be no difference in the cells tagged by changing this option.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/trajectories.jl#L12-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.full_cell-Tuple{Any, Any}' href='#Jutul.full_cell-Tuple{Any, Any}'><span class="jlbinding">Jutul.full_cell</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Inner cell to local cell (full set)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/dd/trivial_map.jl#L35" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_1d_interpolator-Tuple{Any, Any}' href='#Jutul.get_1d_interpolator-Tuple{Any, Any}'><span class="jlbinding">Jutul.get_1d_interpolator</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_1d_interpolator(xs, ys; <keyword arguments>)
```


Get a 1D interpolator `F(x) ≈ y` for a table `xs, ys` that by default does constant extrapolation

**Arguments**
- `xs`: sorted list of parameter points.
  
- `ys`: list of function values with equal length to `xs`
  
- `method=LinearInterpolant`: constructor for the interpolation. Defaults to `LinearInterpolant` which does simple linear interpolation.
  
- `cap_endpoints = true`: Add values so that the endpoints are capped (constant extrapolation). Otherwise, the extrapolation will match the method.
  
- `cap_start = cap_endpoints`: Fine-grained version of cap_endpoints for the start of the interval only (extrapolation for `x < xs[1]`)
  
- `cap_end = cap_endpoints`:Fine-grained version of cap_endpoints for the end of the interval only (extrapolation for `x > xs[end]`)
  

Additional keyword arguments are passed onto the interpolator constructor.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/interpolation.jl#L101-L116" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_2d_interpolator-Tuple{Any, Any, Any}' href='#Jutul.get_2d_interpolator-Tuple{Any, Any, Any}'><span class="jlbinding">Jutul.get_2d_interpolator</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_2d_interpolator(xs, ys, fs; method = BilinearInterpolant, cap_endpoints = true)
```


For `xs` of length `nx` and `ys` of length `ny` generate a 2D interpolation for values given as a `nx` by `ny` matrix. By default `cap_endpoints=true`, and constant extrapolation is used. Fine-grined control over extrapolation can be achieved by setting the keywords arguments `cap_x = (cap_low_x, cap_high_x)` and analogously for `cap_y`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/interpolation.jl#L212-L220" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_ad_entity_scalar-Union{Tuple{T}, Tuple{T, Any}, Tuple{T, Any, Any}} where T<:Real' href='#Jutul.get_ad_entity_scalar-Union{Tuple{T}, Tuple{T, Any}, Tuple{T, Any, Any}} where T<:Real'><span class="jlbinding">Jutul.get_ad_entity_scalar</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_ad_entity_scalar(v::Real, npartials, diag_pos = nothing; <keyword_arguments>)
```


Get scalar with partial derivatives as AD instance.

**Arguments**
- `v::Real`: Value of AD variable.
  
- `npartials`: Number of partial derivatives each AD instance holds.
  
- `diag_pos` = nothing: Position(s) of where to set 1 as the partial derivative instead of zero.
  

**Keyword arguments**
- `tag = nothing`: Tag for AD instance. Two AD values of the different tag cannot interoperate to avoid perturbation confusion (see ForwardDiff documentation).
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L356-L368" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_dependencies-Tuple{Any, Any}' href='#Jutul.get_dependencies-Tuple{Any, Any}'><span class="jlbinding">Jutul.get_dependencies</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Get dependencies of variable when viewed as a secondary variable. Normally autogenerated with @jutul_secondary


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variable_evaluation.jl#L199-L201" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_diagonal_entries-Tuple{JutulEquation, Any}' href='#Jutul.get_diagonal_entries-Tuple{JutulEquation, Any}'><span class="jlbinding">Jutul.get_diagonal_entries</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_diagonal_entries(eq::JutulEquation, eq_s)
```


Get the diagonal entries of a cache, i.e. the entries where entity type and index equals that of the governing equation.

Note: Be very careful about modifications to this array, as this is a view into the internal AD buffers and it is very easy to create inconsistent Jacobians.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L609-L616" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_entity_tag-Tuple{Any, Any}' href='#Jutul.get_entity_tag-Tuple{Any, Any}'><span class="jlbinding">Jutul.get_entity_tag</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_entity_tag(basetag, entity)
```


Combine a base tag (which can be nothing) with a entity to get a tag that captures base tag + entity tag for use with AD initialization.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L454-L459" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_entries-Tuple{CompactAutoDiffCache}' href='#Jutul.get_entries-Tuple{CompactAutoDiffCache}'><span class="jlbinding">Jutul.get_entries</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Get entries of autodiff cache. Entries are AD vectors that hold values and derivatives.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/compact.jl#L1-L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_entries-Tuple{JutulEquation}' href='#Jutul.get_entries-Tuple{JutulEquation}'><span class="jlbinding">Jutul.get_entries</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Get the entries of the main autodiff cache for an equation.

Note: This only gets the .equation field&#39;s entries.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L19-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_mesh_entity_tag-Tuple{JutulMesh, Vararg{Any}}' href='#Jutul.get_mesh_entity_tag-Tuple{JutulMesh, Vararg{Any}}'><span class="jlbinding">Jutul.get_mesh_entity_tag</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_mesh_entity_tag(met::JutulMesh, entity::JutulEntity, tag_group::Symbol, tag_value = missing; throw = true)
```


Get the indices tagged for `entity` in group `tag_group`, optionally for the specific `tag_value`. If `ismissing(tag_value)`, the Dict containing the tag group will be returned.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L1268-L1274" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_parameters-Tuple{SimulationModel}' href='#Jutul.get_parameters-Tuple{SimulationModel}'><span class="jlbinding">Jutul.get_parameters</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_parameters(model::SimulationModel)
```


Get the parameter definitions (as `OrderedDict`) for a given `model`.

Parameters are defined as static values in a forward simulation that combine with the primary variables to compute secondary variables and model equations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L42-L49" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_primary_variable_ordered_entities-Tuple{SimulationModel}' href='#Jutul.get_primary_variable_ordered_entities-Tuple{SimulationModel}'><span class="jlbinding">Jutul.get_primary_variable_ordered_entities</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_primary_variable_ordered_entities(model::SimulationModel)
```


Get only the entities where primary variables are present, sorted by their order in the primary variables.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L186-L190" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_primary_variables-Tuple{SimulationModel}' href='#Jutul.get_primary_variables-Tuple{SimulationModel}'><span class="jlbinding">Jutul.get_primary_variables</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_primary_variables(model::SimulationModel)
```


Get the primary variable definitions (as `OrderedDict`) for a given `model`.

Primary variables are sometimes referred to as solution variables or primary unknowns. The set of primary variables completely determines the state of the system together with the `parameters`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L5-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_secondary_variables-Tuple{SimulationModel}' href='#Jutul.get_secondary_variables-Tuple{SimulationModel}'><span class="jlbinding">Jutul.get_secondary_variables</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_secondary_variables(model::SimulationModel)
```


Get the secondary variable definitions (as `OrderedDict`) for a given `model`.

Secondary variables are variables that can be computed from the primary variables together with the parameters.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L18-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_sparse_arguments-Tuple{Any, Any}' href='#Jutul.get_sparse_arguments-Tuple{Any, Any}'><span class="jlbinding">Jutul.get_sparse_arguments</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_sparse_arguments(storage, model)
```


Get the [`SparsePattern`]@ref for the Jacobian matrix of a given simulator storage and corresponding model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L554-L558" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_tstr' href='#Jutul.get_tstr'><span class="jlbinding">Jutul.get_tstr</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_tstr(dT, lim = 3)
```


Get formatted time string of `dT` given in seconds, limited to `lim` number of units.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/utils.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_variable-Tuple{SimulationModel, Symbol}' href='#Jutul.get_variable-Tuple{SimulationModel, Symbol}'><span class="jlbinding">Jutul.get_variable</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_variable(model::SimulationModel, name::Symbol)
```


Get implementation of variable or parameter with name `name` for the model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L69-L73" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.get_variables-Tuple{SimulationModel}' href='#Jutul.get_variables-Tuple{SimulationModel}'><span class="jlbinding">Jutul.get_variables</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_variables(model::SimulationModel)
```


Get all variable definitions (as `OrderedDict`) for a given `model`.

This is the union of [`get_secondary_variables`](/ref/jutul#Jutul.get_secondary_variables-Tuple{SimulationModel}) and [`get_primary_variables`](/ref/jutul#Jutul.get_primary_variables-Tuple{SimulationModel}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L31-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.global_cell-Tuple{Any, Any}' href='#Jutul.global_cell-Tuple{Any, Any}'><span class="jlbinding">Jutul.global_cell</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Local cell -&gt; global cell (full set)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/dd/trivial_map.jl#L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.global_face-Tuple{Any, Any}' href='#Jutul.global_face-Tuple{Any, Any}'><span class="jlbinding">Jutul.global_face</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Local face -&gt; global face (full set)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/dd/trivial_map.jl#L21" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.initialize_context!-NTuple{4, Any}' href='#Jutul.initialize_context!-NTuple{4, Any}'><span class="jlbinding">Jutul.initialize_context!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Initialize context when setting up a model


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/context.jl#L62-L64" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.initialize_extra_state_fields!-Tuple{Any, JutulModel}' href='#Jutul.initialize_extra_state_fields!-Tuple{Any, JutulModel}'><span class="jlbinding">Jutul.initialize_extra_state_fields!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
initialize_extra_state_fields!(state, model::JutulModel)
```


Add model-dependent changing variables that need to be in state, but are never AD variables themselves (for example status flags).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L269-L274" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.initialize_storage!-Tuple{Any, JutulModel}' href='#Jutul.initialize_storage!-Tuple{Any, JutulModel}'><span class="jlbinding">Jutul.initialize_storage!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
initialize_storage!(storage, model::JutulModel; initialize_state0 = true)
```


Initialize the already allocated storage at the beginning of a simulation. Use this to e.g. set up extra stuff in state0 needed for initializing the simulation loop.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L390-L395" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.initialize_variable_value-Tuple{Any, VectorVariables, AbstractVector}' href='#Jutul.initialize_variable_value-Tuple{Any, VectorVariables, AbstractVector}'><span class="jlbinding">Jutul.initialize_variable_value</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Initializer for the value of non-scalar primary variables


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L335-L337" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.interior_cell-Tuple{Any, Any}' href='#Jutul.interior_cell-Tuple{Any, Any}'><span class="jlbinding">Jutul.interior_cell</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Local cell in full set -&gt; inner cell (or zero)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/dd/trivial_map.jl#L32" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.interpolation_constant_lookup' href='#Jutul.interpolation_constant_lookup'><span class="jlbinding">Jutul.interpolation_constant_lookup</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
interpolation_constant_lookup(X, constant_dx = missing)
```


Generate a lookup table for linear interpolation when dx is evenly spaced.

Note: Setting `constant_dx=true` can lead to incorrect interpolations if the data is not evenly spaced.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/interpolation.jl#L43-L50" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.jutul_message' href='#Jutul.jutul_message'><span class="jlbinding">Jutul.jutul_message</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
jutul_message("Jutul", "Hello, world!", color = :light_blue)
```


Print a line with a colored prefix. The prefix is colored with the `color` argument. The `fancy` argument controls whether the output is colored or not.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/print.jl#L189-L194" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.jutul_output_path' href='#Jutul.jutul_output_path'><span class="jlbinding">Jutul.jutul_output_path</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
pth = jutul_output_path(name = missing; subfolder = "jutul", basedir = missing, create = true)
```


Get path for output. The final path will be found in /basedir/&lt;subfolder/name. If `subfolder=missing`, the path will be set to /basedir/name instead. `name` will be autogenerated if not provided.

Pass the optional input `create = false` to avoid making the directory. To globally set the default output dir, set `ENV["JUTUL_OUTPUT_PATH"]``to your desired`basedir``.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/utils.jl#L939-L948" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.linear_timestep_selection' href='#Jutul.linear_timestep_selection'><span class="jlbinding">Jutul.linear_timestep_selection</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
linear_timestep_selection(x, x0, x1, dt0, dt1)
```


Produce linear estimate of timestep `dt` for some value `x` from observed observations. If the observations have the same `x` or `dt` values, a simple scaling based on the `x1` value is used.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/timesteps.jl#L226-L232" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.load_balanced_endpoint-Tuple{Any, Any, Any}' href='#Jutul.load_balanced_endpoint-Tuple{Any, Any, Any}'><span class="jlbinding">Jutul.load_balanced_endpoint</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
load_balanced_endpoint(block_index, nvals, nblocks)
```


Endpoint for interval `block_index` that subdivides `nvals` into `nblocks` in a load balanced manner. This is done by adding one element to the first set of blocks whenever possible.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/partitioning.jl#L311-L317" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.load_balanced_interval-Tuple{Any, Any, Any}' href='#Jutul.load_balanced_interval-Tuple{Any, Any, Any}'><span class="jlbinding">Jutul.load_balanced_interval</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
load_balanced_interval(b, n, m)
```


Create UnitRange for block b ∈ [1, m] for interval of total length n


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/partitioning.jl#L328-L332" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.local_ad-Tuple{Any, Any, Any}' href='#Jutul.local_ad-Tuple{Any, Any, Any}'><span class="jlbinding">Jutul.local_ad</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
local_ad(state::T, index::I, ad_tag::∂T) where {T, I<:Integer, ∂T}
```


Create local_ad for state for index I of AD tag of type ad_tag


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/local_ad.jl#L180-L184" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.local_cell-Tuple{Any, Any}' href='#Jutul.local_cell-Tuple{Any, Any}'><span class="jlbinding">Jutul.local_cell</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Global cell -&gt; local cell (full set)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/dd/trivial_map.jl#L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.local_face-Tuple{Any, Any}' href='#Jutul.local_face-Tuple{Any, Any}'><span class="jlbinding">Jutul.local_face</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Global face -&gt; local face (full set)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/dd/trivial_map.jl#L30" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.local_residual_view-NTuple{4, Any}' href='#Jutul.local_residual_view-NTuple{4, Any}'><span class="jlbinding">Jutul.local_residual_view</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
local_residual_view(r_buf, model, eq, equation_offset)
```


Get a matrix view of the residual so that, independent of ordering, the column index corresponds to the entity index for the given equation `eq` starting at `equation_offset` in the global residual buffer `r_buf`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L753-L759" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.maximum_value-Tuple{JutulVariables}' href='#Jutul.maximum_value-Tuple{JutulVariables}'><span class="jlbinding">Jutul.maximum_value</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Upper (inclusive) limit for variable.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L98-L100" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.merge_step_report_errors-Tuple{Any}' href='#Jutul.merge_step_report_errors-Tuple{Any}'><span class="jlbinding">Jutul.merge_step_report_errors</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
merge_step_report_errors(data; fn = max)
```


Merge step reports errors of the same type using a pair wise reduction (default: max)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/utils.jl#L1009-L1013" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.mesh_from_gmsh' href='#Jutul.mesh_from_gmsh'><span class="jlbinding">Jutul.mesh_from_gmsh</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
G = mesh_from_gmsh(pth)
G = mesh_from_gmsh()
G = mesh_from_gmsh(pth; verbose = true)
```


Parse a Gmsh file and return a Jutul `UnstructuredMesh` (in 3D only). Requires the Gmsh.jl package to be loaded. If no path is provided in `pth` it is assumed that you are managing the Gmsh state manually and it will use the current selected mesh inside Gmsh. Please note that Gmsh is GPL licensed unless you have obtained another type of license from the authors.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ext/gmsh_ext.jl#L3-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.minimum_value-Tuple{JutulVariables}' href='#Jutul.minimum_value-Tuple{JutulVariables}'><span class="jlbinding">Jutul.minimum_value</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Lower (inclusive) limit for variable.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L103-L105" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.model_accumulation!' href='#Jutul.model_accumulation!'><span class="jlbinding">Jutul.model_accumulation!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
model_accumulation!(acc, sim::HelperSimulator, x, dt = 1.0;
    forces = setup_forces(sim.model),
    update_secondary = true,
    kwarg...
)
```


Compute the accumulation term into Vector acc.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/helper.jl#L106-L114" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.model_residual!' href='#Jutul.model_residual!'><span class="jlbinding">Jutul.model_residual!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
model_residual!(r, sim, x, x0 = missing, dt = 1.0;
    forces = setup_forces(sim.model),
    include_accumulation = true,
    kwarg...
)
```


Fill in the model residual into Vector r.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/helper.jl#L60-L68" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.model_residual-Union{Tuple{T}, Tuple{HelperSimulator{<:Any, <:Any, <:Any, T}, Any, Vararg{Any}}} where T' href='#Jutul.model_residual-Union{Tuple{T}, Tuple{HelperSimulator{<:Any, <:Any, <:Any, T}, Any, Vararg{Any}}} where T'><span class="jlbinding">Jutul.model_residual</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
model_residual(sim::HelperSimulator, x, y = missing; kwarg...)
```


Out of place version of `model_residual!`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/helper.jl#L47-L51" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_boundary_faces' href='#Jutul.number_of_boundary_faces'><span class="jlbinding">Jutul.number_of_boundary_faces</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
number_of_boundary_faces(g)
```


Get the number of boundary/exterior faces in a mesh.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/meshes.jl#L157-L161" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_cells-Tuple{JutulMesh}' href='#Jutul.number_of_cells-Tuple{JutulMesh}'><span class="jlbinding">Jutul.number_of_cells</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
number_of_cells(g)::Integer
```


Get the number of cells in a mesh.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/meshes.jl#L143-L147" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_cells-Tuple{Union{DataDomain, DiscretizedDomain}}' href='#Jutul.number_of_cells-Tuple{Union{DataDomain, DiscretizedDomain}}'><span class="jlbinding">Jutul.number_of_cells</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
number_of_cells(D::Union{DataDomain, DiscretizedDomain})
```


Get the number of cells in a `DataDomain` or `DiscretizedDomain`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/domains.jl#L60-L64" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_degrees_of_freedom-Tuple{JutulModel}' href='#Jutul.number_of_degrees_of_freedom-Tuple{JutulModel}'><span class="jlbinding">Jutul.number_of_degrees_of_freedom</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Total number of degrees of freedom for a model, over all primary variables and all entities.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L14-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_entities-Tuple{Any, JutulEquation}' href='#Jutul.number_of_entities-Tuple{Any, JutulEquation}'><span class="jlbinding">Jutul.number_of_entities</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Get the number of entities (e.g. the number of cells) that the equation is defined on.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L279-L281" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_entities-Tuple{Any, JutulVariables}' href='#Jutul.number_of_entities-Tuple{Any, JutulVariables}'><span class="jlbinding">Jutul.number_of_entities</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Number of entities (e.g. Cells, Faces) a variable is defined on. By default, each primary variable exists on all cells of a discretized domain


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L2-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_entities-Tuple{JutulAutoDiffCache}' href='#Jutul.number_of_entities-Tuple{JutulAutoDiffCache}'><span class="jlbinding">Jutul.number_of_entities</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Get number of entities a cache is defined on.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L3-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_entities-Tuple{T} where T<:AbstractVector' href='#Jutul.number_of_entities-Tuple{T} where T<:AbstractVector'><span class="jlbinding">Jutul.number_of_entities</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Number of entities for vector stored in state (just the number of elements)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L9-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_entities-Tuple{T} where T<:AbstractArray' href='#Jutul.number_of_entities-Tuple{T} where T<:AbstractArray'><span class="jlbinding">Jutul.number_of_entities</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Number of entities for matrix stored in state (convention is number of columns)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L14-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_equations-Tuple{Any, JutulEquation}' href='#Jutul.number_of_equations-Tuple{Any, JutulEquation}'><span class="jlbinding">Jutul.number_of_equations</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Get the total number of equations on the domain of model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L286-L288" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_equations_per_entity-Tuple{JutulModel, JutulEquation}' href='#Jutul.number_of_equations_per_entity-Tuple{JutulModel, JutulEquation}'><span class="jlbinding">Jutul.number_of_equations_per_entity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
n = number_of_equations_per_entity(model::JutulModel, eq::JutulEquation)
```


Get the number of equations per entity. For example, mass balance of two components will have two equations per grid cell (= entity)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L260-L265" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_faces-Tuple{Any}' href='#Jutul.number_of_faces-Tuple{Any}'><span class="jlbinding">Jutul.number_of_faces</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
number_of_faces(g)::Integer
```


Get the number of faces in a mesh.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/meshes.jl#L150-L154" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_faces-Tuple{Union{DataDomain, DiscretizedDomain}}' href='#Jutul.number_of_faces-Tuple{Union{DataDomain, DiscretizedDomain}}'><span class="jlbinding">Jutul.number_of_faces</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
number_of_faces(D::Union{DataDomain, DiscretizedDomain})
```


Get the number of faces in a `DataDomain` or `DiscretizedDomain`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/domains.jl#L69-L73" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_half_faces-Tuple{Union{DataDomain, DiscretizedDomain}}' href='#Jutul.number_of_half_faces-Tuple{Union{DataDomain, DiscretizedDomain}}'><span class="jlbinding">Jutul.number_of_half_faces</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
number_of_half_faces(D::Union{DataDomain, DiscretizedDomain})
```


Get the number of half-faces in a `DataDomain` or `DiscretizedDomain`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/domains.jl#L78-L82" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_partials_per_entity-Tuple{SimulationModel, JutulEntity}' href='#Jutul.number_of_partials_per_entity-Tuple{SimulationModel, JutulEntity}'><span class="jlbinding">Jutul.number_of_partials_per_entity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
number_of_partials_per_entity(model::SimulationModel, entity::JutulEntity)
```


Get the number of local partial derivatives per entity in a `model` for a given [`JutulEntity`](/ref/jutul#Jutul.JutulEntity). This is the sum of [`degrees_of_freedom_per_entity`](/ref/jutul#Jutul.degrees_of_freedom_per_entity-Tuple{Any,%20ScalarVariable}) for all primary variables defined on `entity`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L209-L214" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.number_of_values' href='#Jutul.number_of_values'><span class="jlbinding">Jutul.number_of_values</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



Total number of values for a model, for a given type of variables over all entities


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L34-L36" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.numerical_eltype-Union{Tuple{AbstractArray{T}}, Tuple{T}} where T' href='#Jutul.numerical_eltype-Union{Tuple{AbstractArray{T}}, Tuple{T}} where T'><span class="jlbinding">Jutul.numerical_eltype</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
numerical_eltype(x::AbstractArray{T}) where T
```


Get the numerical eltype (i.e. the inner type of the element type that could potentially be AD)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/local_ad.jl#L112-L117" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.numerical_type-Tuple{T} where T<:Real' href='#Jutul.numerical_type-Tuple{T} where T<:Real'><span class="jlbinding">Jutul.numerical_type</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
numerical_type(::T) where T
```


Get the numerical eltype (i.e. the inner type of the element type that could potentially be AD). This function should be overloaded if you have a custom type that wraps a numeric/potentially AD type.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/local_ad.jl#L126-L132" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.parameters_jacobian_wrt_data_domain-Tuple{Any}' href='#Jutul.parameters_jacobian_wrt_data_domain-Tuple{Any}'><span class="jlbinding">Jutul.parameters_jacobian_wrt_data_domain</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
parameters_jacobian_wrt_data_domain(model; copy = true, config = nothing)
```


Compute the (sparse) Jacobian of parameters with respect to data_domain values (i.e. floating point values). Optionally, `config` can be passed to allow `vectorize_variables` to only include a subset of the parameters.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/vectorization.jl#L204-L210" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.partition' href='#Jutul.partition'><span class="jlbinding">Jutul.partition</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
partition(N::AbstractMatrix, num_coarse, weights = ones(size(N, 2)); partitioner = MetisPartitioner(), groups = nothing, n = maximum(N), group_by_weights = false, buffer_group = true)
```


Partition based on neighborship (with optional groups kept contigious after partitioning)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/partitioning.jl#L239-L244" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.partition_hypergraph' href='#Jutul.partition_hypergraph'><span class="jlbinding">Jutul.partition_hypergraph</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
partition_hypergraph(g, n::Int, partitioner = MetisPartitioner(); expand = true)
```


Partition a hypergraph from [`setup_partitioner_hypergraph`](/ref/jutul#Jutul.setup_partitioner_hypergraph-Tuple{Matrix{Int64}}) using a given partitioner. If the optional `expand` parameter is set to true the result will be expanded to the full graph (i.e. where groups are not condensed).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/partitioning.jl#L433-L439" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.physical_representation-Tuple{Any}' href='#Jutul.physical_representation-Tuple{Any}'><span class="jlbinding">Jutul.physical_representation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
physical_representation(x)
```


Get the physical representation of an object. The physical representation is usually some kind of mesh or domain that represents a physical domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/domains.jl#L7-L12" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.physical_representation-Tuple{DataDomain}' href='#Jutul.physical_representation-Tuple{DataDomain}'><span class="jlbinding">Jutul.physical_representation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
physical_representation(x::DataDomain)
```


Get the underlying physical representation (domain or mesh) that is wrapped.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/domains.jl#L71-L75" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.physical_representation-Tuple{DiscretizedDomain}' href='#Jutul.physical_representation-Tuple{DiscretizedDomain}'><span class="jlbinding">Jutul.physical_representation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
physical_representation(x::DiscretizedDomain)
```


Get the underlying physical representation (domain or mesh) that was discretized.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/domains.jl#L23-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.physical_representation-Tuple{SimulationModel}' href='#Jutul.physical_representation-Tuple{SimulationModel}'><span class="jlbinding">Jutul.physical_representation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
physical_representation(m::SimulationModel)
```


Get the underlying physical representation for the model (domain or mesh)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/core_types/core_types.jl#L319-L323" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.pick_next_timestep-Tuple{IterationTimestepSelector, Vararg{Any, 9}}' href='#Jutul.pick_next_timestep-Tuple{IterationTimestepSelector, Vararg{Any, 9}}'><span class="jlbinding">Jutul.pick_next_timestep</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
pick_next_timestep(sel::IterationTimestepSelector, sim, config, dt_prev, dT, forces, reports, current_reports, step_index, new_step)
```


Pick the next time-step for `IterationTimestepSelector`. This function uses the number of iterations from previous timesteps to estimate the relationship between the last and the new time step.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/timesteps.jl#L61-L67" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.plot_cell_data!-Tuple' href='#Jutul.plot_cell_data!-Tuple'><span class="jlbinding">Jutul.plot_cell_data!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_cell_data!(ax, mesh, data; kwarg...)
```


Mutating version of `plot_cell_data` that plots into an existing Makie `Axis`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ext/makie_ext.jl#L123-L127" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.plot_cell_data-Tuple' href='#Jutul.plot_cell_data-Tuple'><span class="jlbinding">Jutul.plot_cell_data</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_cell_data(mesh::JutulMesh, data::Vector; kwarg...)
plot_cell_data(mesh, data;
    cells = nothing,
    faces = nothing,
    boundaryfaces = nothing
)
```


Plot cell-wise values (as a vector) on the mesh. Optionally, indices `cells`, `faces` or `boundaryfaces` can be passed to limit the plotting to a specific selection of entities.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ext/makie_ext.jl#L102-L113" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.plot_interactive-Tuple' href='#Jutul.plot_interactive-Tuple'><span class="jlbinding">Jutul.plot_interactive</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_interactive(mesh, vector_of_dicts; kwarg...)
```


Launch an interactive plot of a mesh with the given `vector_of_dicts` (or just a dict). Each dict can have cell data either as vectors (one value per cell) or matrices (one column per cell).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ext/makie_ext.jl#L8-L14" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.plot_mesh!-Tuple' href='#Jutul.plot_mesh!-Tuple'><span class="jlbinding">Jutul.plot_mesh!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_mesh!(ax, mesh)
```


Mutating version of `plot_mesh` that plots into an existing Makie `Axis` instance.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ext/makie_ext.jl#L57-L62" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.plot_mesh-Tuple' href='#Jutul.plot_mesh-Tuple'><span class="jlbinding">Jutul.plot_mesh</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_mesh(mesh)
plot_mesh(mesh;
    cells = nothing,
    faces = nothing,
    boundaryfaces = nothing,
    outer = false,
    color = :lightblue,
)
```


Plot a `mesh` with uniform colors. Optionally, indices `cells`, `faces` or `boundaryfaces` can be passed to limit the plotting to a specific selection of entities.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ext/makie_ext.jl#L33-L46" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.plot_mesh_edges!-Tuple' href='#Jutul.plot_mesh_edges!-Tuple'><span class="jlbinding">Jutul.plot_mesh_edges!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_mesh_edges!(ax, mesh; kwarg...)
```


Plot the edges of all cells on the exterior of a mesh into existing Makie `Axis` `ax`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ext/makie_ext.jl#L87-L92" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.plot_mesh_edges-Tuple' href='#Jutul.plot_mesh_edges-Tuple'><span class="jlbinding">Jutul.plot_mesh_edges</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_mesh_edges(mesh; kwarg...)
```


Plot the edges of all cells on the exterior of a mesh.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ext/makie_ext.jl#L73-L77" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.prepare_step_storage-Tuple{Any, Any, Missing}' href='#Jutul.prepare_step_storage-Tuple{Any, Any, Missing}'><span class="jlbinding">Jutul.prepare_step_storage</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
prepare_step_storage(storage, model, ::Missing)
```


Initialize storage for prepare_step_handler.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/simulator.jl#L64-L68" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.process_partition' href='#Jutul.process_partition'><span class="jlbinding">Jutul.process_partition</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
p = process_partition(g::JutulMesh, partition; weights = missing)
p = process_partition(N::Vector{Tuple{Int, Int}}, partition; weights = missing)
p = process_partition(N::Matrix{Int}, partition; weights = missing)
```


Perform processing of `partition` on mesh or neighborship to make sure that all coarse blocks in the partition are connected. Optionally, weights can be passed that will be used to remove connections with weights that are zero (or within floating point precision of zero).

The resulting partition vector will be a copy and will have additional coarse blocks inserted when blocks were split into two or more components.

```julia
g = CartesianMesh((5,))
p = [1, 1, 2, 1, 1] # Block 1 is disconnected!
process_partition(g, p)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/partitioning.jl#L101-L119" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.read_results-Tuple{Any}' href='#Jutul.read_results-Tuple{Any}'><span class="jlbinding">Jutul.read_results</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



states, reports = read_results(pth; read_states = true, read_reports = true)

Read results from a given `output_path` provded to `simulate` or `simulator_config`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/utils.jl#L608-L612" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.relative_increment_limit-Tuple{JutulVariables}' href='#Jutul.relative_increment_limit-Tuple{JutulVariables}'><span class="jlbinding">Jutul.relative_increment_limit</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Relative allowable change for variable during a nonlinear update. A variable with value |x| and relative limit 0.2 cannot change more than |x|*0.2.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L91-L95" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.replace_variables!-Tuple{Any}' href='#Jutul.replace_variables!-Tuple{Any}'><span class="jlbinding">Jutul.replace_variables!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
replace_variables!(model, throw = true, varname = vardef, varname2 = vardef2)
```


Replace one or more variables that already exists in the model (primary, secondary or parameters) with a new definition.

**Arguments**
- `model`: instance where variables is to be replaced
  
- `varname=vardef::JutulVariables`: replace variable with `varname` by `vardef`
  
- `throw=true`: throw an error if the named variable definition is not found in primary or secondary, otherwise silently return the `model`.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L137-L148" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.scalarize_primary_variable-Tuple{Any, Any, ScalarVariable, Any}' href='#Jutul.scalarize_primary_variable-Tuple{Any, Any, ScalarVariable, Any}'><span class="jlbinding">Jutul.scalarize_primary_variable</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scalarize_primary_variable(model, source_vec, var::Jutul.ScalarVariable, index)
```


Scalarize a primary variable. For scalars, this means getting the value itself.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/scalarization.jl#L40-L44" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.scalarize_primary_variables' href='#Jutul.scalarize_primary_variables'><span class="jlbinding">Jutul.scalarize_primary_variables</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
scalarize_primary_variables(model, state, pvars = model.primary_variables)
```


Create a vector where each entry corresponds to a tuple of values that minimally defines the given variables. All variables must belong to the same type of entity. This is checked by this function.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/scalarization.jl#L91-L97" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.scalarize_primary_variables!-Union{Tuple{T}, Tuple{Array{Jutul.ScalarizedJutulVariables{T}, 1}, Any, Any, NamedTuple}} where T' href='#Jutul.scalarize_primary_variables!-Union{Tuple{T}, Tuple{Array{Jutul.ScalarizedJutulVariables{T}, 1}, Any, Any, NamedTuple}} where T'><span class="jlbinding">Jutul.scalarize_primary_variables!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scalarize_primary_variables!(V::Vector{T}, model, state, pvars::NamedTuple) where T
```


Scalarize into array. See [`scalarize_primary_variables`](/ref/jutul#Jutul.scalarize_primary_variables) for more details.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/scalarization.jl#L115-L119" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.scalarized_primary_variable_type-Tuple{Any, ScalarVariable}' href='#Jutul.scalarized_primary_variable_type-Tuple{Any, ScalarVariable}'><span class="jlbinding">Jutul.scalarized_primary_variable_type</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scalarized_primary_variable_type(model, var::Jutul.ScalarVariable)
```


Get the type of a scalarized numerical variable (=Float64 for variables that are already represented as scalars)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/scalarization.jl#L17-L22" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.set_default_tolerances-Tuple{Any}' href='#Jutul.set_default_tolerances-Tuple{Any}'><span class="jlbinding">Jutul.set_default_tolerances</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
set_default_tolerances(model)
```


Set default tolerances for the nonlinear convergence check of the governing equations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L771-L775" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.set_parameters!-Tuple{Any}' href='#Jutul.set_parameters!-Tuple{Any}'><span class="jlbinding">Jutul.set_parameters!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
set_parameters!(model, parname = pardef)
```


Set a parameter with name `varname` to the definition `vardef` (adding if it does not exist)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L119-L123" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.set_primary_variables!-Tuple{Any}' href='#Jutul.set_primary_variables!-Tuple{Any}'><span class="jlbinding">Jutul.set_primary_variables!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
set_primary_variables!(model, varname = vardef)
set_primary_variables!(model, varname1 = vardef1, varname2 = vardef2)
```


Set a primary variable with name `varname` to the definition `vardef` (adding if it does not exist)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L96-L101" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.set_secondary_variables!-Tuple{Any}' href='#Jutul.set_secondary_variables!-Tuple{Any}'><span class="jlbinding">Jutul.set_secondary_variables!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
set_secondary_variables!(model, varname = vardef)
set_secondary_variables!(model, varname1 = vardef1, varname2 = vardef2)
```


Set a secondary variable with name `varname` to the definition `vardef` (adding if it does not exist)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L107-L113" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.setup_adjoint_storage-Tuple{Any}' href='#Jutul.setup_adjoint_storage-Tuple{Any}'><span class="jlbinding">Jutul.setup_adjoint_storage</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setup_adjoint_storage(model; state0 = setup_state(model), parameters = setup_parameters(model))
```


Set up storage for use with `solve_adjoint_sensitivities!`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/gradients.jl#L103-L107" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.setup_forces-Tuple{JutulModel}' href='#Jutul.setup_forces-Tuple{JutulModel}'><span class="jlbinding">Jutul.setup_forces</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setup_forces(model::JutulModel; force_name = force_value)
```


Set up forces for a given model. Keyword arguments varies depending on what the model supports.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L876-L881" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.setup_parameter_optimization-Tuple{Any, Any, Any, Any, Any, Any, Vararg{Any}}' href='#Jutul.setup_parameter_optimization-Tuple{Any, Any, Any, Any, Any, Any, Vararg{Any}}'><span class="jlbinding">Jutul.setup_parameter_optimization</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setup_parameter_optimization(model, state0, param, dt, forces, G, opt_cfg = optimization_config(model, param);
                                                        grad_type = :adjoint,
                                                        config = nothing,
                                                        print = 1,
                                                        copy_case = true,
                                                        param_obj = false,
                                                        kwarg...)
```


Set up function handles for optimizing the case defined by the inputs to `simulate` together with a per-timestep objective function `G`. The objective should be on the form of sum over all steps, where each element of the sum is evaluated by `model, state, dt, step_no, forces`.

Generally calling either of the functions will mutate the data Dict. The options are: F_o(x) -&gt; evaluate objective dF_o(dFdx, x) -&gt; evaluate gradient of objective, mutating dFdx (may trigger evaluation of F_o) F_and_dF(F, dFdx, x) -&gt; evaluate F and/or dF. Value of nothing will mean that the corresponding entry is skipped.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/optimization.jl#L17-L36" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.setup_parameters-Tuple{DataDomain, JutulModel}' href='#Jutul.setup_parameters-Tuple{DataDomain, JutulModel}'><span class="jlbinding">Jutul.setup_parameters</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setup_parameters(model::JutulModel; name = value)
```


Set up a parameter storage for a given model with values for the parameter defined in the model.

**Arguments**
- `name=value`: The name of the parameter together with the value(s) of the parameter.
  

A scalar (or short vector of the right size for [`VectorVariables`](/ref/jutul#Jutul.VectorVariables)) will be repeated over the entire domain, while a vector (or matrix for [`VectorVariables`](/ref/jutul#Jutul.VectorVariables)) with length (number of columns for [`VectorVariables`](/ref/jutul#Jutul.VectorVariables)) equal to the entity count (for example, number of cells for a cell variable) will be used directly.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L293-L303" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.setup_partitioner_hypergraph-Tuple{Matrix{Int64}}' href='#Jutul.setup_partitioner_hypergraph-Tuple{Matrix{Int64}}'><span class="jlbinding">Jutul.setup_partitioner_hypergraph</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setup_partitioner_hypergraph(N::Matrix{Int};
    num_nodes::Int = maximum(N),
    num_edges::Int = size(N, 2),
    node_weights::Vector{Int} = ones(Int, num_nodes),
    edge_weights::Vector{Int} = ones(Int, num_edges),
    groups = [Int[]]
)
```


Set up a hypergraph structure for a given neighborship matrix. `N` should be a matrix with two rows, with one pair of cells in each column. Optionally node and edge weights can be provided. If a list of groups are provided, these nodes will be accumulated together in the hypergraph.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/partitioning.jl#L339-L352" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.setup_state!' href='#Jutul.setup_state!'><span class="jlbinding">Jutul.setup_state!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
setup_state!(state, model::JutulModel, init_values::AbstractDict = Dict())
```


Initialize primary variables and other state fields, given initial values as a Dict


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L254-L258" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.setup_state-Tuple{JutulModel, Vararg{Any}}' href='#Jutul.setup_state-Tuple{JutulModel, Vararg{Any}}'><span class="jlbinding">Jutul.setup_state</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setup_state(model::JutulModel, name1 = value1, name2 = value2)
```


Set up a state for a given model with values for the primary variables defined in the model. Normally all primary variables must be initialized in this way.

**Arguments**
- `name=value`: The name of the primary variable together with the value(s) used to initialize the primary variable.
  

A scalar (or short vector of the right size for [`VectorVariables`](/ref/jutul#Jutul.VectorVariables)) will be repeated over the entire domain, while a vector (or matrix for [`VectorVariables`](/ref/jutul#Jutul.VectorVariables)) with length (number of columns for [`VectorVariables`](/ref/jutul#Jutul.VectorVariables)) equal to the entity count (for example, number of cells for a cell variable) will be used directly.

Note: You likely want to overload [`setup_state!`]@ref for a custom model instead of `setup_state`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L226-L239" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.setup_state_and_parameters-Tuple{DataDomain, JutulModel, AbstractDict}' href='#Jutul.setup_state_and_parameters-Tuple{DataDomain, JutulModel, AbstractDict}'><span class="jlbinding">Jutul.setup_state_and_parameters</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
state, prm = setup_state_and_parameters(model, init)
```


Simultaneously set up state and parameters from a single `init` file (typically a `Dict` containing values that might either be initial values or parameters)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L327-L332" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.setup_storage!-Tuple{Any, JutulModel}' href='#Jutul.setup_storage!-Tuple{Any, JutulModel}'><span class="jlbinding">Jutul.setup_storage!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setup_storage!(storage, model::JutulModel; setup_linearized_system = true,
                                                setup_equations = true,
                                                state0 = setup_state(model),
                                                parameters = setup_parameters(model),
                                                tag = nothing,
                                                state0_ad = false,
                                                state_ad = true,
                                                kwarg...)
```


Allocate storage for a given model. The storage consists of all dynamic quantities used in the simulation. The default implementation allocates properties, equations and linearized system.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L407-L419" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.setup_storage-Tuple{JutulModel}' href='#Jutul.setup_storage-Tuple{JutulModel}'><span class="jlbinding">Jutul.setup_storage</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setup_storage(model::JutulModel; kwarg...)
```


Allocate storage for the model. You should overload setup_storage! if you have a custom definition.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L378-L383" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.si_unit-Tuple{Symbol}' href='#Jutul.si_unit-Tuple{Symbol}'><span class="jlbinding">Jutul.si_unit</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
si_unit(u::Union{String, Symbol})
```


Get the multiplicative SI unit conversion factor for a single unit. The return value is given so that `x*si_unit(:name)` will convert `x` to the SI representation of the unit with the given name.

**Examples**

```julia
julia> si_unit(:day) # Get days represented as seconds
86400.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/units/interface.jl#L75-L87" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.si_units-Tuple' href='#Jutul.si_units-Tuple'><span class="jlbinding">Jutul.si_units</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
u1_val = si_units(u1)
meter = si_units(:meter)
meter, hour = si_units(:meter, :hour)
```


Get multiplicative SI unit conversion factors for multiple units simultaneously. The return value will be a `Tuple` of values, one for each input argument. Each input arguments can be either a `String` a `Symbol`.

**Examples**

```julia
julia> meter, hour = si_units(:meter, :hour)
(1.0, 3600.0)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/units/interface.jl#L56-L70" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.simulate!-Tuple{JutulSimulator, AbstractVector}' href='#Jutul.simulate!-Tuple{JutulSimulator, AbstractVector}'><span class="jlbinding">Jutul.simulate!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
simulate!(sim::JutulSimulator, timesteps::AbstractVector;
    forces = nothing,
    config = nothing,
    initialize = true,
    restart = nothing,
    state0 = nothing,
    parameters = nothing,
    kwarg...
)
```


Non-allocating (or perhaps less allocating) version of [`simulate!`](/ref/jutul#Jutul.simulate!-Tuple{JutulSimulator,%20AbstractVector}).

**Arguments**
- `initialize=true`: Perform internal updates as if this is the first time 
  

See also [`simulate`](/ref/jutul#Jutul.simulate-Tuple{Any,%20JutulModel,%20AbstractVector}) for additional supported input arguments.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/simulator.jl#L131-L148" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.simulate-Tuple{Any, JutulModel, AbstractVector}' href='#Jutul.simulate-Tuple{Any, JutulModel, AbstractVector}'><span class="jlbinding">Jutul.simulate</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
simulate(state0, model, timesteps, parameters = setup_parameters(model))
simulate(state0, model, timesteps, info_level = 3)
simulate(state0, model, timesteps; <keyword arguments>)
```


Simulate a set of `timesteps` with `model` for the given initial `state0` and optionally specific parameters. Additional keyword arguments are passed onto [`simulator_config`](/ref/jutul#Jutul.simulator_config-Tuple{Any}) and [`simulate!`](/ref/jutul#Jutul.simulate!-Tuple{JutulSimulator,%20AbstractVector}). This interface is primarily for convenience, as all storage for the simulator is allocated upon use and discared upon return. If you want to perform multiple simulations with the same model it is advised to instead instantiate [`Simulator`](/ref/jutul#Jutul.Simulator-Tuple{Any})  and combine it with [`simulate!`](/ref/jutul#Jutul.simulate!-Tuple{JutulSimulator,%20AbstractVector}).

**Arguments**
- `state0::Dict`: initial state, typically created using [`setup_state`](/ref/jutul#Jutul.setup_state-Tuple{JutulModel,%20Vararg{Any}}) for the `model` in use.
  
- `model::JutulModel`: model that describes the discretized system to solve, for example [`SimulationModel`](/ref/jutul#Jutul.SimulationModel-Tuple{Any,%20Any}) or [`MultiModel`](/ref/jutul#Jutul.MultiModel).
  
- `timesteps::AbstractVector`: Vector of desired report steps. The simulator will perform time integration until `sum(timesteps)`  is reached, providing outputs at the end of each report step.
  
- `parameters=setup_parameters(model)`: Optional overrides the default parameters for the model.
  
- `forces=nothing`: Either `nothing` (for no forces), a single set of forces from `setup_forces(model)` or a `Vector` of such forces with equal length to `timesteps`.
  
- `restart=nothing`: If an integer is provided, the simulation will attempt to restart from that step. Requires that `output_path` is provided here or in the `config`.
  
- `config=simulator_config(model)`: Configuration `Dict` that holds many fine grained settings for output, linear solver, time-steps, outputs etc.
  

Additional arguments are passed onto [`simulator_config`](/ref/jutul#Jutul.simulator_config-Tuple{Any}).

See also [`simulate!`](/ref/jutul#Jutul.simulate!-Tuple{JutulSimulator,%20AbstractVector}), [`Simulator`](/ref/jutul#Jutul.Simulator-Tuple{Any}), [`SimulationModel`](/ref/jutul#Jutul.SimulationModel-Tuple{Any,%20Any}), [`simulator_config`](/ref/jutul#Jutul.simulator_config-Tuple{Any}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/simulator.jl#L84-L109" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.simulate-Tuple{Any, JutulSimulator, AbstractVector}' href='#Jutul.simulate-Tuple{Any, JutulSimulator, AbstractVector}'><span class="jlbinding">Jutul.simulate</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
simulate(state0, sim::JutulSimulator, timesteps::AbstractVector; parameters = nothing, kwarg...)
```


Simulate a set of `timesteps` with `simulator` for the given initial `state0` and optionally specific parameters.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/simulator.jl#L120-L124" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.simulator_config-Tuple{Any}' href='#Jutul.simulator_config-Tuple{Any}'><span class="jlbinding">Jutul.simulator_config</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
simulator_config(sim; info_level = 3, linear_solver = GenericKrylov())
```


Set up a simulator configuration object that can be passed onto [`simulate!`](/ref/jutul#Jutul.simulate!-Tuple{JutulSimulator,%20AbstractVector}).

There are many options available to configure a given simulator. The best way to get an overview of these possible configuration options is to instatiate the config without any arguments and inspecting the resulting table by calling `simulator_config(sim)` in the REPL.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/config.jl#L96-L106" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.solve_adjoint_forces-Tuple{JutulCase, Any, Any}' href='#Jutul.solve_adjoint_forces-Tuple{JutulCase, Any, Any}'><span class="jlbinding">Jutul.solve_adjoint_forces</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
solve_adjoint_forces(case::JutulCase, res::SimResult, G)
```


Solve the adjoint equations for the forces in `case` given the simulation result `res` and the objective function `G`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/force_gradients.jl#L395-L400" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.solve_adjoint_sensitivities!-NTuple{6, Any}' href='#Jutul.solve_adjoint_sensitivities!-NTuple{6, Any}'><span class="jlbinding">Jutul.solve_adjoint_sensitivities!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
solve_adjoint_sensitivities!(∇G, storage, states, state0, timesteps, G; forces = setup_forces(model))
```


Non-allocating version of `solve_adjoint_sensitivities`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/gradients.jl#L222-L226" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.solve_adjoint_sensitivities-NTuple{4, Any}' href='#Jutul.solve_adjoint_sensitivities-NTuple{4, Any}'><span class="jlbinding">Jutul.solve_adjoint_sensitivities</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
solve_adjoint_sensitivities(model, states, reports_or_timesteps, G; extra_timing = false, state0 = setup_state(model), forces = setup_forces(model), raw_output = false, kwarg...)
```


Compute sensitivities of `model` parameter with name `target` for objective function `G`.

The objective function is at the moment assumed to be a sum over all states on the form: `obj = Σₙ G(model, state, dt_n, n, forces_for_step_n)`

Solves the adjoint equations: For model equations F the gradient with respect to parameters p is     ∇ₚG = Σₙ (∂Fₙ / ∂p)ᵀ λₙ where n ∈ [1, N]. Given Lagrange multipliers λₙ from the adjoint equations     (∂Fₙ / ∂xₙ)ᵀ λₙ = - (∂J / ∂xₙ)ᵀ - (∂Fₙ₊₁ / ∂xₙ)ᵀ λₙ₊₁ where the last term is omitted for step n = N and G is the objective function.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/gradients.jl#L3-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.solve_numerical_sensitivities-NTuple{5, Any}' href='#Jutul.solve_numerical_sensitivities-NTuple{5, Any}'><span class="jlbinding">Jutul.solve_numerical_sensitivities</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
solve_numerical_sensitivities(model, states, reports, G, target;
                                            forces = setup_forces(model),
                                            state0 = setup_state(model),
                                            parameters = setup_parameters(model),
                                            epsilon = 1e-8)
```


Compute sensitivities of `model` parameter with name `target` for objective function `G`.

This method uses numerical perturbation and is primarily intended for testing of `solve_adjoint_sensitivities`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/gradients.jl#L606-L616" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.solve_timestep!-NTuple{5, Any}' href='#Jutul.solve_timestep!-NTuple{5, Any}'><span class="jlbinding">Jutul.solve_timestep!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
solve_timestep!(sim, dT, forces, max_its, config; <keyword arguments>)
```


Internal function for solving a single time-step with fixed driving forces.

**Arguments**
- `sim`: `Simulator` instance.
  
- `dT`: time-step to be solved
  
- `forces`: Driving forces for the time-step
  
- `max_its`: Maximum number of steps/Newton iterations.
  
- `config`: Configuration for solver (typically output from `simulator_config`).
  

Note: This function is exported for fine-grained simulation workflows. The general [`simulate`](/ref/jutul#Jutul.simulate-Tuple{Any,%20JutulModel,%20AbstractVector}) interface is both easier to use and performs additional validation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/simulator/simulator.jl#L253-L267" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.subvariable-Tuple{Any, Any}' href='#Jutul.subvariable-Tuple{Any, Any}'><span class="jlbinding">Jutul.subvariable</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
subvariable(var, map)
```


Get subvariable of Jutul variable


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/dd/submodels.jl#L144-L148" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.successful_reports' href='#Jutul.successful_reports'><span class="jlbinding">Jutul.successful_reports</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
successful_reports(reports; step = length(reports), n = 1)
```


Get last n successful reports starting at the end of `step` and reversing backwards until `n` values have been found. `n` can be set to `Inf` to produce all successful reports.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/timesteps.jl#L211-L217" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.successful_reports-2' href='#Jutul.successful_reports-2'><span class="jlbinding">Jutul.successful_reports</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
successful_reports(old_reports, current_reports, step_index, n = 1)
```


Get the `n` last successful solve reports from all previous reports (old_reports) and the current ministep set. You can optionally provide a function that replaces the default definition of `success=r->r[:success]`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/timesteps.jl#L174-L180" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.synchronize-Tuple{JutulContext}' href='#Jutul.synchronize-Tuple{JutulContext}'><span class="jlbinding">Jutul.synchronize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Synchronize backend after allocations.

Some backends may require notification that storage has been allocated.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/context.jl#L66-L71" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.tpfv_geometry' href='#Jutul.tpfv_geometry'><span class="jlbinding">Jutul.tpfv_geometry</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
tpfv_geometry(g)
```


Generate two-point finite-volume geometry for a given grid, if supported.

See also [`TwoPointFiniteVolumeGeometry`](/ref/jutul#Jutul.TwoPointFiniteVolumeGeometry).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/meshes.jl#L10-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.transfer-Tuple{Any, Any}' href='#Jutul.transfer-Tuple{Any, Any}'><span class="jlbinding">Jutul.transfer</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Transfer v to the representation expected by a given context.

For the defalt context, the transfer function does nothing. For other context such as the CUDA version, it may convert integers and floats to other types (e.g. Float32) and Arrays to CuArrays.

You will likely have to implement some transfer operators for your own types if you want to simulate with a non-default context.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/context.jl#L1-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.two_point_potential_drop-NTuple{5, Real}' href='#Jutul.two_point_potential_drop-NTuple{5, Real}'><span class="jlbinding">Jutul.two_point_potential_drop</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Two-point potential drop with gravity (generic)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/conservation/flux.jl#L297-L299" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.unsafe_reinterpret-Tuple{Any, Any, Any}' href='#Jutul.unsafe_reinterpret-Tuple{Any, Any, Any}'><span class="jlbinding">Jutul.unsafe_reinterpret</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
unsafe_reinterpret(Vt, v, n)
```


Unsafely reinterpret v as a n length vector of value type Vt


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/linsolve/utils.jl#L95-L99" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_before_step!-NTuple{4, Any}' href='#Jutul.update_before_step!-NTuple{4, Any}'><span class="jlbinding">Jutul.update_before_step!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>




<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L931-L933" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_equation!-Tuple{Any, JutulEquation, Any, Any, Any}' href='#Jutul.update_equation!-Tuple{Any, JutulEquation, Any, Any, Any}'><span class="jlbinding">Jutul.update_equation!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Update equation based on currently stored properties


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L512-L514" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_equations!' href='#Jutul.update_equations!'><span class="jlbinding">Jutul.update_equations!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
update_equations!(storage, model, dt = nothing)
```


Update the governing equations using the current set of primary variables, parameters and secondary variables. Does not fill linearized system.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L710-L714" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_equations_and_apply_forces!-NTuple{4, Any}' href='#Jutul.update_equations_and_apply_forces!-NTuple{4, Any}'><span class="jlbinding">Jutul.update_equations_and_apply_forces!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
update_equations_and_apply_forces!(storage, model, dt, forces; time = NaN)
```


Update the model equations and apply boundary conditions and forces. Does not fill linearized system.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L692-L696" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_linearized_system!' href='#Jutul.update_linearized_system!'><span class="jlbinding">Jutul.update_linearized_system!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
update_linearized_system!(storage, model::JutulModel; <keyword arguments>)
```


Update the linearized system with the current set of equations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L725-L729" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_linearized_system_equation!-Tuple{AbstractArray, Any, Any, JutulEquation, CompactAutoDiffCache}' href='#Jutul.update_linearized_system_equation!-Tuple{AbstractArray, Any, Any, JutulEquation, CompactAutoDiffCache}'><span class="jlbinding">Jutul.update_linearized_system_equation!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Update a linearized system based on the values and derivatives in the equation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/equations.jl#L490-L492" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_parameter_before_step!-NTuple{6, Any}' href='#Jutul.update_parameter_before_step!-NTuple{6, Any}'><span class="jlbinding">Jutul.update_parameter_before_step!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
update_parameter_before_step!(prm_val, prm, storage, model, dt, forces)
```


Update parameters before time-step. Used for hysteretic parameters.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L977-L981" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_secondary_variable!' href='#Jutul.update_secondary_variable!'><span class="jlbinding">Jutul.update_secondary_variable!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



Update a secondary variable. Normally autogenerated with @jutul_secondary


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variable_evaluation.jl#L207-L209" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_state_dependents!-Tuple{Any, JutulModel, Any, Any}' href='#Jutul.update_state_dependents!-Tuple{Any, JutulModel, Any, Any}'><span class="jlbinding">Jutul.update_state_dependents!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
update_state_dependents!(storage, model, dt, forces; time = NaN, update_secondary = true)
```


Perform updates of everything that depends on the state: A full linearization for the current primary variables.

This includes properties, governing equations and the linearized system itself.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/models.jl#L677-L683" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_values!-Tuple{AbstractArray, AbstractArray}' href='#Jutul.update_values!-Tuple{AbstractArray, AbstractArray}'><span class="jlbinding">Jutul.update_values!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
update_values!(x, dx)
```


Replace values (for non-Real types, direct assignment)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L393-L397" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.update_values!-Tuple{AbstractArray{<:ForwardDiff.Dual}, AbstractArray{<:Real}}' href='#Jutul.update_values!-Tuple{AbstractArray{<:ForwardDiff.Dual}, AbstractArray{<:Real}}'><span class="jlbinding">Jutul.update_values!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
update_values!(x, dx)
```


Replace values of `x` in-place by `y`, leaving `x` with the values of y and the partials of `x`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L379-L383" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.value-Tuple{AbstractDict}' href='#Jutul.value-Tuple{AbstractDict}'><span class="jlbinding">Jutul.value</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
value(d::Dict)
```


Call value on all elements of some Dict.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L425-L428" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.value-Tuple{Any}' href='#Jutul.value-Tuple{Any}'><span class="jlbinding">Jutul.value</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Take value of AD.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ad/ad.jl#L414-L416" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.values_per_entity-Tuple{Any, JutulVariables}' href='#Jutul.values_per_entity-Tuple{Any, JutulVariables}'><span class="jlbinding">Jutul.values_per_entity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Number of values held by a primary variable. Normally this is equal to the number of degrees of freedom, but some special primary variables are most conveniently defined by having N values and N-1 independent variables.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L79-L82" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.variable_scale-Tuple{JutulVariables}' href='#Jutul.variable_scale-Tuple{JutulVariables}'><span class="jlbinding">Jutul.variable_scale</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Define a &quot;typical&quot; numerical value for a variable to scale the linear system entries.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variables/utils.jl#L190-L192" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.write_reports_to_mat_format' href='#Jutul.write_reports_to_mat_format'><span class="jlbinding">Jutul.write_reports_to_mat_format</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
write_reports_to_mat_format(reports::Vector, path = jutul_output_path(); name = "report", config = missing, verbose = false)
```


Write the reports to MAT files named &quot;report_1&quot;, &quot;report_2&quot;, ... to the given path.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/utils.jl#L1098-L1102" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.@jutul_secondary-Tuple{Any}' href='#Jutul.@jutul_secondary-Tuple{Any}'><span class="jlbinding">Jutul.@jutul_secondary</span></a> <Badge type="info" class="jlObjectType jlMacro" text="Macro" /></summary>



Designate the function as updating a secondary variable.

A generic evaluator is then defined, together with a function for getting the dependencies of that function upon the state. This is most easily documented with an example. If we define the following function annotated with the macro when updating the array containing the values of `MyVarType` realized for some model:

```julia
@jutul_secondary function some_fn!(target, var::MyVarType, model, a, b, c, ix)
    for i in ix
        target[i] = a[i] + b[i] / c[i]
    end
end
```


The purpose of the macro is to translate this into two functions. The first defines for the dependencies of the function with respect to the fields of the state (primary variables, secondary variables and parameters):

```julia
function get_dependencies(var::MyVarType, model)
   return (:a, :b, :c)
end
```


The second function defines a generic version that takes in state, and automatically expands the set of dependencies into `getfield` calls.

```julia
function update_secondary_variable!(array_target, var::MyVarType, model, state, ix)
    some_fn!(array_target, var, model, state.a, state.b, state.c, ix)
end
```


Note that the input names of arguments 4 to end-1 matter, as these will be fetched from state, exactly as written.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/variable_evaluation.jl#L4-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.ConvergenceMonitorRelaxation-Tuple{}' href='#Jutul.ConvergenceMonitors.ConvergenceMonitorRelaxation-Tuple{}'><span class="jlbinding">Jutul.ConvergenceMonitors.ConvergenceMonitorRelaxation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ConvergenceMonitorRelaxation(; w_min = 0.1, dw = 0.2, dw_increase = nothing, dw_decrease = nothing, w_max = 1.0)
```


Relaxation strategy based on convergence monitoring. Requires that the  simulation has been configured with a `ConvergenceMonitorCuttingCriterion` so that the convergence status is available in the reports. See corresponding iplementation of `select_nonlinear_relaxation_model` below for details.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/relaxation.jl#L8-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.compute_contraction_factor-Tuple{Any, Any}' href='#Jutul.ConvergenceMonitors.compute_contraction_factor-Tuple{Any, Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.compute_contraction_factor</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
compute_contraction_factor(r, N)
```


Compute contraction factor from a number of Newton iterate distances from convergence (defined with a user-prescribed distance metric). For more than two iterates, the contraction factor is estimated using least-squares assuming the iterates follow a geometric series: rₙ = Θⁿr₀. The function also computes the target contraction factor under this assumption assuming N iterations left.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/contraction_factors.jl#L1-L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.compute_distance-Tuple{Any}' href='#Jutul.ConvergenceMonitors.compute_distance-Tuple{Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.compute_distance</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
compute_distance(report; distance_function = r -> scaled_residual_norm(r), mapping = v -> maximum(v))
```


Compute distance from convergence using a user-defined distance function, and optionally apply a mapping to the distance. The function returns the distance and the names of equation residual norms used in the distance computation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/distance_functions.jl#L1-L7" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.flatten_dict' href='#Jutul.ConvergenceMonitors.flatten_dict'><span class="jlbinding">Jutul.ConvergenceMonitors.flatten_dict</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
flatten_dict(input_dict::Dict, separator::String = ".", trail = [])
```


Flatten a dict of dicts into a vector of values and a vector of names. The names are on the format `"key1<separator>key2<separator>key3"` and the values are the corresponding values in the dict.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/utils.jl#L87-L93" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.get_model_residuals-Tuple{Any}' href='#Jutul.ConvergenceMonitors.get_model_residuals-Tuple{Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.get_model_residuals</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_model_residuals(report)
```


Get all residual norms for all equations of a model from a nonlinear iteration report, scaled by their respective toelrances. The function returns a dict of dicts with the residual norms for each equation on the form

```
residuals[equation_type][equation][norm] = res/tol
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/utils.jl#L27-L36" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.get_multimodel_residuals-Tuple{Any}' href='#Jutul.ConvergenceMonitors.get_multimodel_residuals-Tuple{Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.get_multimodel_residuals</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_multimodel_residuals(report)
```


Get all residual norms for all equations of all models from a nonlinear iteration report, scaled by their respective toelrances. The function returns a dict of dicts with the residual norms for each equation on the form

```
residuals[model][equation_type][equation][norm] = res/tol
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/utils.jl#L3-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.make_report-NTuple{4, Any}' href='#Jutul.ConvergenceMonitors.make_report-NTuple{4, Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.make_report</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
make_report(θ, θ_target, oscillation, status)
```


Utility for generating a report from the convergence monitor.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/cutting_criterions.jl#L147-L151" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.nonconverged_equations-Tuple{Any}' href='#Jutul.ConvergenceMonitors.nonconverged_equations-Tuple{Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.nonconverged_equations</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
nonconverged_equations(report)
```


Compute distance to convergence as non-converged equations, e.g., 1.0 for non-converged and 0.0 for converged. The final distance is typically taken as the sum of the output from this function.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/distance_functions.jl#L38-L44" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.oscillation' href='#Jutul.ConvergenceMonitors.oscillation'><span class="jlbinding">Jutul.ConvergenceMonitors.oscillation</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
oscillation(contraction_factors, tol)
```


Check if the contraction factors are oscillating. The function checks if the contraction factors are oscillating by checking if the last three contraction factors are below a user-defined tolerance defining &quot;slow&quot; convergece. The function returns true if the contraction factors are oscillating, and false otherwise. If less than three contraction factors are provided, the function returns false.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/contraction_factors.jl#L27-L36" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.print_convergence_status-Tuple{Jutul.ConvergenceMonitors.ConvergenceMonitorCuttingCriterion, Any, Any}' href='#Jutul.ConvergenceMonitors.print_convergence_status-Tuple{Jutul.ConvergenceMonitors.ConvergenceMonitorCuttingCriterion, Any, Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.print_convergence_status</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
print_convergence_status(cc::ConvergenceMonitorCuttingCriterion, it, it0)
```


Utility for printing the status of the convergence monitor.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/cutting_criterions.jl#L163-L167" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.process_name-Tuple{Any}' href='#Jutul.ConvergenceMonitors.process_name-Tuple{Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.process_name</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
process_name(name)
```


Process a names to be suibale as dictionary keys.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/utils.jl#L74-L78" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.reset!-Tuple{Jutul.ConvergenceMonitors.ConvergenceMonitorCuttingCriterion, Any, Any}' href='#Jutul.ConvergenceMonitors.reset!-Tuple{Jutul.ConvergenceMonitors.ConvergenceMonitorCuttingCriterion, Any, Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.reset!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
reset!(cc::ConvergenceMonitorCuttingCriterion, template, max_iter)
```


Utility for resetting convergence monitor cutting criterion.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/cutting_criterions.jl#L120-L124" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.scaled_residual_norm-Tuple{Any}' href='#Jutul.ConvergenceMonitors.scaled_residual_norm-Tuple{Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.scaled_residual_norm</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scaled_residual_norm(report)
```


Compute distance to convergence as the residual norm of the equations scaled by their respective tolerances.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/distance_functions.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.set_convergence_monitor_cutting_criterion!-Tuple{Any}' href='#Jutul.ConvergenceMonitors.set_convergence_monitor_cutting_criterion!-Tuple{Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.set_convergence_monitor_cutting_criterion!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
set_convergence_monitor_cutting_criterion!(config; max_nonlinear_iterations = 50, kwargs...)
```


Utility for setting `ConvergenceMonitorCuttingCriterion` to the simulator config. The function also adjusts the maximum number of nonlinear iterations (to 50 by default).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/cutting_criterions.jl#L20-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.ConvergenceMonitors.set_convergence_monitor_relaxation!-Tuple{Any}' href='#Jutul.ConvergenceMonitors.set_convergence_monitor_relaxation!-Tuple{Any}'><span class="jlbinding">Jutul.ConvergenceMonitors.set_convergence_monitor_relaxation!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
set_convergence_monitor_relaxation!(config; max_nonlinear_iterations = 50, convergence_monitor_args = NamedTuple(), relaxation_args...)
```


Utility for setting `ConvergenceMonitorRelaxation` to the simulator config. This function also sets the a `ConvergenceMonitorCuttingCriterion` to the config – see `set_convergence_monitor_cutting_criterion!`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/relaxation.jl#L26-L32" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.cutting_criterion-Tuple{Jutul.ConvergenceMonitors.ConvergenceMonitorCuttingCriterion, Vararg{Any, 9}}' href='#Jutul.cutting_criterion-Tuple{Jutul.ConvergenceMonitors.ConvergenceMonitorCuttingCriterion, Vararg{Any, 9}}'><span class="jlbinding">Jutul.cutting_criterion</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
Jutul.cutting_criterion(cc::ConvergenceMonitorCuttingCriterion, sim, dt, forces, it, max_iter, cfg, e, step_reports, relaxation)
```


Cutting ctriterion based on monitoring convergence. The function computes the contraction factor from the distance to convergence (in a user-defined metric) between consecutive iterates. The function also computes the target contraction factor under the assumption that the iterates follow a geometric series, and checks if the iterates are oscillating. Based on this, the iterate is classified as &quot;good&quot;, &quot;ok&quot;, or &quot;bad&quot;, and a counter for number of violations is updated accordingly (+1 for &quot;bad&quot;, -1 for &quot;good&quot;). The timestep is aborted if the number of violations exceeds a user-defined limit.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/cutting_criterions.jl#L43-L54" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Jutul.select_nonlinear_relaxation_model-Tuple{Any, Jutul.ConvergenceMonitors.ConvergenceMonitorRelaxation, Any, Any}' href='#Jutul.select_nonlinear_relaxation_model-Tuple{Any, Jutul.ConvergenceMonitors.ConvergenceMonitorRelaxation, Any, Any}'><span class="jlbinding">Jutul.select_nonlinear_relaxation_model</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
Jutul.select_nonlinear_relaxation_model(model, rel_type::ConvergenceMonitorRelaxation, reports, ω)
```


Relaxation strategy based on convergence monitoring. The relaxation factor is decreased if the convergence status is &quot;bad&quot;, and increased if the status is &quot;good&quot;, determined by contraction factors computed from distance to convergence (in a user-defined metric) between consecutive iterates. These are computed during the cutting criterion check and stored in the reports.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/ConvergenceMonitors/relaxation.jl#L48-L56" target="_blank" rel="noreferrer">source</a></Badge>

</details>

