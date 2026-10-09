abstract type AbstractEmbeddedFracture end

const Point3 = SVector{3, Float64}
const Point2 = SVector{2, Float64}

"""Physical properties of one fracture, in SI units."""
struct FractureProperties
    aperture::Float64
    permeability::Float64
    porosity::Float64
end

function FractureProperties(;
        aperture = 5e-4,
        hydraulic_aperture = aperture,
        permeability = hydraulic_aperture^2/12,
        porosity = 0.5
    )
    isfinite(aperture) && aperture > 0 || throw(ArgumentError("aperture must be positive and finite"))
    isfinite(hydraulic_aperture) && hydraulic_aperture >= 0 || throw(ArgumentError("hydraulic_aperture must be non-negative and finite"))
    isfinite(permeability) && permeability >= 0 || throw(ArgumentError("permeability must be non-negative and finite"))
    isfinite(porosity) && 0 < porosity <= 1 || throw(ArgumentError("porosity must lie in (0, 1]"))
    return FractureProperties(aperture, permeability, porosity)
end

"""
    PolygonFracture(vertices; aperture=5e-4, permeability=aperture^2/12, porosity=0.5)

A simple planar polygon in 3D. Vertices are supplied in boundary order, as a
vector of points or a 3-by-n matrix. Concave polygons are supported; holes and
self-intersections are rejected. A repeated closing vertex is optional.
"""
struct PolygonFracture <: AbstractEmbeddedFracture
    vertices::Vector{Point3}
    normal::Point3
    properties::FractureProperties
end

function PolygonFracture(vertices; kwarg...)
    if vertices isa AbstractMatrix
        size(vertices, 1) == 3 || throw(ArgumentError("vertices must be a 3-by-n matrix"))
        pts = [Point3(vertices[:, i]) for i in axes(vertices, 2)]
    else
        pts = Point3.(vertices)
    end
    if length(pts) > 1 && first(pts) == last(pts)
        pop!(pts)
    end
    length(pts) >= 3 || throw(ArgumentError("a polygon requires at least three vertices"))
    all(p -> all(isfinite, p), pts) || throw(ArgumentError("vertices must be finite"))
    origin = first(pts)
    n = zero(Point3)
    for i in 2:length(pts)-1
        n += cross(pts[i] - origin, pts[i+1] - origin)
    end
    scale = maximum(p -> norm(p - origin), pts)
    norm(n) > 1e-12*scale^2 || throw(ArgumentError("polygon has zero area"))
    n /= norm(n)
    tol = max(1e-10*scale, 64eps(maximum(p -> maximum(abs, p), pts)))
    all(p -> abs(dot(n, p - origin)) <= tol, pts) || throw(ArgumentError("polygon is not planar"))
    u, v = plane_basis(n)
    pts2 = [Point2(dot(p-origin, u), dot(p-origin, v)) for p in pts]
    validate_polygon(pts2, tol)
    return PolygonFracture(pts, n, FractureProperties(; kwarg...))
end

"""
    RectangleFracture(center, normal, lengths; direction=nothing, angle=0.0, kwarg...)

A planar rectangle with full side lengths `lengths`. `direction` specifies the
first in-plane axis (projected onto the plane); `angle` rotates the axes in
radians about `normal`. Physical properties use the same keywords as
[`PolygonFracture`](@ref).
"""
struct RectangleFracture <: AbstractEmbeddedFracture
    center::Point3
    normal::Point3
    axes::Tuple{Point3, Point3}
    lengths::Tuple{Float64, Float64}
    properties::FractureProperties
end

function RectangleFracture(center, normal, lengths; direction = nothing, angle = 0.0, kwarg...)
    c, n, axes = fracture_frame(center, normal, direction, angle)
    length(lengths) == 2 || throw(ArgumentError("lengths must contain two full side lengths"))
    all(x -> isfinite(x) && x > 0, lengths) || throw(ArgumentError("lengths must be positive and finite"))
    return RectangleFracture(c, n, axes, Tuple(Float64.(lengths)), FractureProperties(; kwarg...))
end

"""
    EllipticalFracture(center, normal, radii; direction=nothing, angle=0.0, segments=64, kwarg...)

A parametrized planar ellipse (a flattened ellipsoid) with two semi-axes
`radii`. Geometry is approximated by an inscribed polygon with `segments`
vertices. Increase `segments` for tighter area accuracy. `EllipsoidFracture`
is an alias for this planar fracture type.
"""
struct EllipticalFracture <: AbstractEmbeddedFracture
    center::Point3
    normal::Point3
    axes::Tuple{Point3, Point3}
    radii::Tuple{Float64, Float64}
    segments::Int
    properties::FractureProperties
end

const EllipsoidFracture = EllipticalFracture

function EllipticalFracture(center, normal, radii; direction = nothing, angle = 0.0, segments = 64, kwarg...)
    c, n, axes = fracture_frame(center, normal, direction, angle)
    length(radii) == 2 || throw(ArgumentError("radii must contain two semi-axis lengths"))
    all(x -> isfinite(x) && x > 0, radii) || throw(ArgumentError("radii must be positive and finite"))
    segments >= 8 || throw(ArgumentError("segments must be at least eight"))
    return EllipticalFracture(c, n, axes, Tuple(Float64.(radii)), segments, FractureProperties(; kwarg...))
end

function plane_basis(n)
    axis = argmin(abs.(n))
    e = setindex(zero(Point3), 1.0, axis)
    u = normalize(e-dot(e, n)*n)
    return u, cross(n, u)
end

function fracture_frame(center, normal, direction, angle)
    c, n = Point3(center), Point3(normal)
    all(isfinite, c) && all(isfinite, n) && norm(n) > 0 || throw(ArgumentError("center and nonzero normal must be finite"))
    isfinite(angle) || throw(ArgumentError("angle must be finite"))
    n = normalize(n)
    if isnothing(direction)
        u, v = plane_basis(n)
    else
        u = Point3(direction)
        u -= dot(u, n)*n
        all(isfinite, u) && norm(u) > 1e-12 || throw(ArgumentError("direction must have a nonzero component in the fracture plane"))
        u = normalize(u)
        v = cross(n, u)
    end
    return c, n, (cos(angle)*u + sin(angle)*v, -sin(angle)*u + cos(angle)*v)
end

fracture_polygon(f::PolygonFracture) = copy(f.vertices)

function fracture_polygon(f::RectangleFracture)
    u, v = f.axes
    a, b = f.lengths ./ 2
    return [f.center - a*u - b*v, f.center + a*u - b*v,
            f.center + a*u + b*v, f.center - a*u + b*v]
end

function fracture_polygon(f::EllipticalFracture)
    u, v = f.axes
    a, b = f.radii
    return [f.center + a*cos(t)*u + b*sin(t)*v for t in range(0, 2pi; length = f.segments+1)[1:end-1]]
end
