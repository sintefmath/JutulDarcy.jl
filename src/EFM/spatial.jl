"""
One node of an axis-aligned bounding-box tree, with coordinates in type `T`.

Fields:
- `low`, `high`: Inclusive 3D bounds in metres for all cells below this node.
- `left`, `right`: Child indices in `CellSpatialIndex.nodes`; both are zero for
  a leaf.
- `first`, `last`: Inclusive range into `CellSpatialIndex.cells` for this node.
"""
struct BoundingBoxNode{T <: AbstractFloat}
    low::Point3{T}
    high::Point3{T}
    left::Int
    right::Int
    first::Int
    last::Int
end

"""
Bounding-box tree reused for matrix clipping, fracture intersections and well
completion searches. Queries return sorted local cell IDs.

Fields:
- `nodes`: Tree nodes; node 1 is the root, or the vector is empty for no cells.
- `cells`: Permutation of local cell IDs defining contiguous tree ranges.
- `low`, `high`: Inclusive bounds in metres, indexed by local cell ID rather
  than by the permutation. Matrix bounds use Jutul's `cell_bounding_box`.
"""
struct CellSpatialIndex{T <: AbstractFloat}
    nodes::Vector{BoundingBoxNode{T}}
    cells::Vector{Int}
    low::Vector{Point3{T}}
    high::Vector{Point3{T}}
end

function CellSpatialIndex(low::Vector{Point3{T}}, high::Vector{Point3{T}}; leaf_size = 16) where T
    cells = collect(eachindex(low))
    nodes = BoundingBoxNode{T}[]
    function build(first, last)
        lo, hi = low[cells[first]], high[cells[first]]
        for i in first+1:last
            lo = min.(lo, low[cells[i]])
            hi = max.(hi, high[cells[i]])
        end
        node = length(nodes)+1
        push!(nodes, BoundingBoxNode(lo, hi, 0, 0, first, last))
        if last-first+1 > leaf_size
            axis = argmax(hi-lo)
            sort!(view(cells, first:last); by = c -> low[c][axis]+high[c][axis])
            mid = (first+last)÷2
            left, right = build(first, mid), build(mid+1, last)
            nodes[node] = BoundingBoxNode(lo, hi, left, right, first, last)
        end
        return node
    end
    if !isempty(cells)
        build(1, length(cells))
    end
    return CellSpatialIndex(nodes, cells, low, high)
end

function CellSpatialIndex(mesh::UnstructuredMesh{3}, ::Type{T} = Jutul.float_type(mesh)) where T
    bounds = [Jutul.CutCellMeshes.cell_bounding_box(mesh, c) for c in 1:number_of_cells(mesh)]
    return CellSpatialIndex(Point3{T}[first(b) for b in bounds], Point3{T}[last(b) for b in bounds])
end

boxes_overlap(al, ah, bl, bh, tol) = all(al .<= bh .+ tol) && all(bl .<= ah .+ tol)

function query_cells(index, low, high, tol; origin = nothing, normal = nothing)
    out = Int[]
    isempty(index.nodes) && return out
    stack = [1]
    while !isempty(stack)
        node = index.nodes[pop!(stack)]
        boxes_overlap(node.low, node.high, low, high, tol) || continue
        if !isnothing(normal)
            center, radius = (node.high+node.low)/2, (node.high-node.low)/2
            abs(dot(center-origin, normal)) <= dot(radius, abs.(normal))+tol || continue
        end
        if node.left == 0
            for i in node.first:node.last
                cell = index.cells[i]
                if boxes_overlap(index.low[cell], index.high[cell], low, high, tol)
                    if !isnothing(normal)
                        center = (index.high[cell]+index.low[cell])/2
                        radius = (index.high[cell]-index.low[cell])/2
                        abs(dot(center-origin, normal)) <= dot(radius, abs.(normal))+tol || continue
                    end
                    push!(out, cell)
                end
            end
        else
            push!(stack, node.left, node.right)
        end
    end
    sort!(out)
    return out
end
