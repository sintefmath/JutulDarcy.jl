"""
Abstract planar fracture with coordinates and physical properties stored in
floating-point type `T`. Concrete types provide `normal`, `properties`, and
`fracture_polygon(fracture)` for both EFM and cut-cell DFM preprocessing.
"""
abstract type AbstractEmbeddedFracture{T <: AbstractFloat} end

const Point3{T} = SVector{3, T}
const Point2{T} = SVector{2, T}

"""
Physical properties of one fracture, stored in floating-point type `T`.

Fields:
- `aperture`: Mechanical aperture in metres, used for physical storage.
- `permeability`: Isotropic in-plane permeability in m². The constructor
  defaults to `hydraulic_aperture^2/12`; hydraulic aperture is not stored.
- `porosity`: Fracture pore-volume fraction in `(0, 1]`.
"""
struct FractureProperties{T <: AbstractFloat}
    aperture::T
    permeability::T
    porosity::T
end

function FractureProperties(;
        float_type = Float64,
        aperture = 5e-4,
        hydraulic_aperture = aperture,
        permeability = hydraulic_aperture^2/12,
        porosity = 0.5
    )
    isfinite(aperture) && aperture > 0 || throw(ArgumentError("aperture must be positive and finite"))
    isfinite(hydraulic_aperture) && hydraulic_aperture >= 0 || throw(ArgumentError("hydraulic_aperture must be non-negative and finite"))
    isfinite(permeability) && permeability >= 0 || throw(ArgumentError("permeability must be non-negative and finite"))
    isfinite(porosity) && 0 < porosity <= 1 || throw(ArgumentError("porosity must lie in (0, 1]"))
    T = float_type
    return FractureProperties{T}(T(aperture), T(permeability), T(porosity))
end

"""
    PolygonFracture(vertices; aperture=5e-4, permeability=aperture^2/12, porosity=0.5)

A simple planar polygon in 3D. Vertices are supplied in boundary order, as a
vector of points or a 3-by-n matrix. Concave polygons are supported; holes and
self-intersections are rejected. A repeated closing vertex is optional.

Fields:
- `vertices`: Boundary-ordered 3D coordinates in metres, without a closing copy.
- `normal`: Unit normal following the vertex winding.
- `properties`: Mechanical aperture, in-plane permeability and porosity.
- `mesh_options`: EFM internal cell counts or target cell lengths. Constructor
  keywords are `resolution` and `cell_length`; omit both for matrix clipping.
"""
struct PolygonFracture{T <: AbstractFloat} <: AbstractEmbeddedFracture{T}
    vertices::Vector{Point3{T}}
    normal::Point3{T}
    properties::FractureProperties{T}
    mesh_options::FractureMeshOptions{T}
end

function PolygonFracture(vertices; float_type = nothing, resolution = nothing, cell_length = nothing, kwarg...)
    if isnothing(float_type)
        if vertices isa AbstractMatrix
            float_type = typeof(float(first(vertices)))
        else
            float_type = typeof(float(first(first(vertices))))
        end
    end
    T = float_type
    if vertices isa AbstractMatrix
        size(vertices, 1) == 3 || throw(ArgumentError("vertices must be a 3-by-n matrix"))
        pts = [Point3{T}(vertices[:, i]) for i in axes(vertices, 2)]
    else
        pts = Point3{T}.(vertices)
    end
    if length(pts) > 1 && first(pts) == last(pts)
        pop!(pts)
    end
    length(pts) >= 3 || throw(ArgumentError("a polygon requires at least three vertices"))
    all(p -> all(isfinite, p), pts) || throw(ArgumentError("vertices must be finite"))
    origin = first(pts)
    n = zero(Point3{T})
    for i in 2:length(pts)-1
        n += cross(pts[i] - origin, pts[i+1] - origin)
    end
    scale = maximum(p -> norm(p - origin), pts)
    norm(n) > 1e-12*scale^2 || throw(ArgumentError("polygon has zero area"))
    n /= norm(n)
    tol = max(1e-10*scale, 64eps(maximum(p -> maximum(abs, p), pts)))
    all(p -> abs(dot(n, p - origin)) <= tol, pts) || throw(ArgumentError("polygon is not planar"))
    u, v = plane_basis(n)
    pts2 = [Point2{T}(dot(p-origin, u), dot(p-origin, v)) for p in pts]
    validate_polygon(pts2, tol)
    options = FractureMeshOptions(; resolution = resolution, cell_length = cell_length, float_type = T)
    return PolygonFracture(pts, n, FractureProperties(; float_type = T, kwarg...), options)
end

"""
    RectangleFracture(center, normal, lengths; direction=nothing, angle=0.0, kwarg...)

A planar rectangle with full side lengths `lengths`. `direction` specifies the
first in-plane axis (projected onto the plane); `angle` rotates the axes in
radians about `normal`. Physical properties use the same keywords as
[`PolygonFracture`](@ref).

Fields:
- `center`: Rectangle centre in 3D, in metres.
- `normal`: Unit plane normal.
- `axes`: Orthonormal in-plane axes `(u, v)`, with `cross(u, v) == normal`.
- `lengths`: Full side lengths along `(u, v)`, in metres.
- `properties`: Mechanical aperture, in-plane permeability and porosity.
- `mesh_options`: EFM internal cell counts or target cell lengths, supplied
  through `resolution` or `cell_length`.
"""
struct RectangleFracture{T <: AbstractFloat} <: AbstractEmbeddedFracture{T}
    center::Point3{T}
    normal::Point3{T}
    axes::Tuple{Point3{T}, Point3{T}}
    lengths::Tuple{T, T}
    properties::FractureProperties{T}
    mesh_options::FractureMeshOptions{T}
end

function RectangleFracture(center, normal, lengths; direction = nothing, angle = 0.0,
        float_type = typeof(float(first(center))), resolution = nothing, cell_length = nothing, kwarg...)
    T = float_type
    c, n, axes = fracture_frame(center, normal, direction, angle, T)
    length(lengths) == 2 || throw(ArgumentError("lengths must contain two full side lengths"))
    all(x -> isfinite(x) && x > 0, lengths) || throw(ArgumentError("lengths must be positive and finite"))
    options = FractureMeshOptions(; resolution = resolution, cell_length = cell_length, float_type = T)
    return RectangleFracture(c, n, axes, Tuple(T.(lengths)), FractureProperties(; float_type = T, kwarg...), options)
end

"""
    EllipticalFracture(center, normal, radii; direction=nothing, angle=0.0, segments=64, kwarg...)

A parametrized planar ellipse (a flattened ellipsoid) with two semi-axes
`radii`. Geometry is approximated by an inscribed polygon with `segments`
vertices. Increase `segments` for tighter area accuracy. `EllipsoidFracture`
is an alias for this planar fracture type.

Fields:
- `center`: Ellipse centre in 3D, in metres.
- `normal`: Unit plane normal.
- `axes`: Orthonormal in-plane axes `(u, v)`, with `cross(u, v) == normal`.
- `radii`: Semi-axis lengths along `(u, v)`, in metres.
- `segments`: Number of vertices in the inscribed boundary polygon.
- `properties`: Mechanical aperture, in-plane permeability and porosity.
- `mesh_options`: EFM internal cell counts or target cell lengths, supplied
  through `resolution` or `cell_length`. Independent of boundary `segments`.
"""
struct EllipticalFracture{T <: AbstractFloat} <: AbstractEmbeddedFracture{T}
    center::Point3{T}
    normal::Point3{T}
    axes::Tuple{Point3{T}, Point3{T}}
    radii::Tuple{T, T}
    segments::Int
    properties::FractureProperties{T}
    mesh_options::FractureMeshOptions{T}
end

const EllipsoidFracture = EllipticalFracture

function EllipticalFracture(center, normal, radii; direction = nothing, angle = 0.0, segments = 64,
        float_type = typeof(float(first(center))), resolution = nothing, cell_length = nothing, kwarg...)
    T = float_type
    c, n, axes = fracture_frame(center, normal, direction, angle, T)
    length(radii) == 2 || throw(ArgumentError("radii must contain two semi-axis lengths"))
    all(x -> isfinite(x) && x > 0, radii) || throw(ArgumentError("radii must be positive and finite"))
    segments >= 8 || throw(ArgumentError("segments must be at least eight"))
    options = FractureMeshOptions(; resolution = resolution, cell_length = cell_length, float_type = T)
    return EllipticalFracture(c, n, axes, Tuple(T.(radii)), segments, FractureProperties(; float_type = T, kwarg...), options)
end

# This basis also fixes the default orientation of rectangular/elliptical axes.
# Keep it stable; cut-cell projection uses a different reference convention.
function plane_basis(n::Point3{T}) where T
    axis = argmin(abs.(n))
    e = setindex(zero(Point3{T}), 1.0, axis)
    u = normalize(e-dot(e, n)*n)
    return u, cross(n, u)
end

function fracture_frame(center, normal, direction, angle, ::Type{T}) where T
    angle = T(angle)
    c, n = Point3{T}(center), Point3{T}(normal)
    all(isfinite, c) && all(isfinite, n) && norm(n) > 0 || throw(ArgumentError("center and nonzero normal must be finite"))
    isfinite(angle) || throw(ArgumentError("angle must be finite"))
    n = normalize(n)
    if isnothing(direction)
        u, v = plane_basis(n)
    else
        u = Point3{T}(direction)
        u -= dot(u, n)*n
        all(isfinite, u) && norm(u) > 1e-12 || throw(ArgumentError("direction must have a nonzero component in the fracture plane"))
        u = normalize(u)
        v = cross(n, u)
    end
    return c, n, (cos(angle)*u + sin(angle)*v, -sin(angle)*u + cos(angle)*v)
end

"""
    fracture_polygon(fracture::AbstractEmbeddedFracture)

Return the fracture boundary as a new vector of 3D points in metres, in the
fracture's floating-point type. Vertices follow the plane normal's winding and
do not repeat the first vertex at the end. Ellipses use their configured
inscribed polygon. This boundary is shared by EFM and cut-cell DFM preprocessing.
"""
fracture_polygon(f::PolygonFracture) = copy(f.vertices)

function fracture_polygon(f::RectangleFracture)
    u, v = f.axes
    a, b = f.lengths ./ 2
    return [f.center - a*u - b*v, f.center + a*u - b*v,
            f.center + a*u + b*v, f.center - a*u + b*v]
end

function fracture_polygon(f::EllipticalFracture{T}) where T
    u, v = f.axes
    a, b = f.radii
    return [f.center + a*cos(t)*u + b*sin(t)*v for t in range(zero(T), 2T(pi); length = f.segments+1)[1:end-1]]
end
