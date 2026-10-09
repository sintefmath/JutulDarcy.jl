struct BoundingBoxNode
    low::Point3
    high::Point3
    left::Int
    right::Int
    first::Int
    last::Int
end

struct CellSpatialIndex
    nodes::Vector{BoundingBoxNode}
    cells::Vector{Int}
    low::Vector{Point3}
    high::Vector{Point3}
end

function CellSpatialIndex(low::Vector{Point3}, high::Vector{Point3}; leaf_size = 16)
    cells = collect(eachindex(low))
    nodes = BoundingBoxNode[]
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

function CellSpatialIndex(mesh::UnstructuredMesh{3})
    bounds = [Jutul.CutCellMeshes.cell_bounding_box(mesh, c) for c in 1:number_of_cells(mesh)]
    return CellSpatialIndex(Point3[first(b) for b in bounds], Point3[last(b) for b in bounds])
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
