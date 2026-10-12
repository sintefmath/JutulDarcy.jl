"""
    FractureMeshOptions(; resolution=nothing, cell_length=nothing, float_type=Float64)

Internal EFM mesh settings for one fracture. Omit both options to use the
matrix-clipped cells. Specify one of `resolution` or `cell_length` to use an
independent planar grid clipped to the fracture and active matrix cells.
`resolution=:matrix` explicitly selects matrix clipping.

Fields:
- `resolution`: `(nu, nv)` positive cell counts along the two fracture axes,
  or `nothing`. An integer constructor argument sets both counts.
- `cell_length`: Positive target grid spacing in metres along both axes,
  stored as a pair in floating-point type `T`, or `nothing`. A scalar sets
  both spacings. Actual spacings are no larger than the requested lengths.

Rectangle and ellipse grids follow their configured axes. Polygon grids use
their `FractureFrame` axes. Resolution counts describe the full bounding grid;
clipping can remove cells or split them into disconnected components.
These settings affect EFM; cut-cell DFM retains its conforming discretization.
"""
struct FractureMeshOptions{T <: AbstractFloat}
    resolution::Union{Nothing, Tuple{Int, Int}}
    cell_length::Union{Nothing, Tuple{T, T}}
end

function FractureMeshOptions(; resolution = nothing, cell_length = nothing, float_type = Float64)
    T = float_type
    T <: AbstractFloat || throw(ArgumentError("float_type must be an AbstractFloat type"))
    isnothing(resolution) || isnothing(cell_length) || throw(ArgumentError("specify resolution or cell_length, not both"))
    if resolution === :matrix
        resolution = nothing
    end
    if !isnothing(resolution)
        if resolution isa Integer
            resolution = (resolution, resolution)
        end
        resolution isa Tuple && length(resolution) == 2 || throw(ArgumentError("resolution must be a positive integer or a pair of positive integers"))
        all(n -> n isa Integer && !(n isa Bool) && n > 0, resolution) || throw(ArgumentError("resolution counts must be positive integers"))
        resolution = (Int(resolution[1]), Int(resolution[2]))
    end
    if !isnothing(cell_length)
        if cell_length isa Number
            cell_length = (cell_length, cell_length)
        end
        cell_length isa Tuple && length(cell_length) == 2 || throw(ArgumentError("cell_length must be a positive length or a pair of positive lengths"))
        all(h -> h isa Real && isfinite(h) && h > 0, cell_length) || throw(ArgumentError("cell_length values must be positive and finite"))
        cell_length = (T(cell_length[1]), T(cell_length[2]))
        all(h -> isfinite(h) && h > 0, cell_length) || throw(ArgumentError("cell_length cannot be represented as a positive finite value in $T"))
    end
    return FractureMeshOptions{T}(resolution, cell_length)
end
