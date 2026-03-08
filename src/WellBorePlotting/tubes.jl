# Contents of this file is mostly AI generated.

# --------- Small vector helpers ----------
@inline norm3(v::NTuple{3,Float64}) = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
@inline function normalize3(v::NTuple{3,Float64})
    n = norm3(v)
    n == 0.0 && return (0.0, 0.0, 0.0)
    (v[1]/n, v[2]/n, v[3]/n)
end
@inline cross(a::NTuple{3,Float64}, b::NTuple{3,Float64}) =
    (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline dot(a::NTuple{3,Float64}, b::NTuple{3,Float64}) =
    a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

# Rotate vector v around unit axis k by angle theta (Rodrigues)
function rotate_around_axis(v::NTuple{3,Float64}, k::NTuple{3,Float64}, θ::Float64)
    c, s = cos(θ), sin(θ)
    kv = dot(k, v)
    cx = cross(k, v)
    (
        v[1]*c + cx[1]*s + k[1]*kv*(1-c),
        v[2]*c + cx[2]*s + k[2]*kv*(1-c),
        v[3]*c + cx[3]*s + k[3]*kv*(1-c),
    )
end

# --------- Catmull-Rom spline sampling (C^1) ----------
"""
    sample_centerline(points; samples_per_seg=8)

Upsamples a polyline `points` (2×M or 3×M) to a smooth Catmull–Rom spline.
Returns `(P::Vector{NTuple{3,Float64}}, T::Vector{NTuple{3,Float64}}, seg_of_sample::Vector{Int})`
with positions `P`, unit tangents `T`, and the originating segment index (1-based) per sample.
"""
function sample_centerline(points::AbstractMatrix{<:Real}; samples_per_seg::Int=8)
    D, M = size(points)
    @assert D == 2 || D == 3 "points must be 2×M or 3×M"
    # Pad to 3D if needed
    P3 = Tuple{Float64,Float64,Float64}[]
    for j in 1:M
        if D == 2
            push!(P3, (float(points[1,j]), float(points[2,j]), 0.0))
        else
            push!(P3, (float(points[1,j]), float(points[2,j]), float(points[3,j])))
        end
    end

    # Tangents for Catmull–Rom (chord-length, interior points)
    function tangent(i)
        if i == 1
            normalize3((P3[2][1]-P3[1][1], P3[2][2]-P3[1][2], P3[2][3]-P3[1][3]))
        elseif i == M
            normalize3((P3[M][1]-P3[M-1][1], P3[M][2]-P3[M-1][2], P3[M][3]-P3[M-1][3]))
        else
            u = (P3[i+1][1]-P3[i-1][1], P3[i+1][2]-P3[i-1][2], P3[i+1][3]-P3[i-1][3])
            normalize3(u)
        end
    end

    P = NTuple{3,Float64}[]
    T = NTuple{3,Float64}[]
    seg_of = Int[]

    # Cubic Hermite basis
    H00(t) = 2t^3 - 3t^2 + 1
    H10(t) = t^3 - 2t^2 + t
    H01(t) = -2t^3 + 3t^2
    H11(t) = t^3 - t^2

    for i in 1:(M-1)
        p0, p1 = P3[i], P3[i+1]
        m0 = tangent(i)
        m1 = tangent(i+1)

        # Scale tangents by segment length for “natural” speed
        seg = (p1[1]-p0[1], p1[2]-p0[2], p1[3]-p0[3])
        L   = max(norm3(seg), eps())
        m0s = (m0[1]*L, m0[2]*L, m0[3]*L)
        m1s = (m1[1]*L, m1[2]*L, m1[3]*L)

        for k in 0:(i==M-1 ? samples_per_seg : samples_per_seg-1)
            t = k/samples_per_seg
            # Position
            px = H00(t)*p0[1] + H10(t)*m0s[1] + H01(t)*p1[1] + H11(t)*m1s[1]
            py = H00(t)*p0[2] + H10(t)*m0s[2] + H01(t)*p1[2] + H11(t)*m1s[2]
            pz = H00(t)*p0[3] + H10(t)*m0s[3] + H01(t)*p1[3] + H11(t)*m1s[3]
            # Tangent (derivative of Hermite)
            dH00 = 6t^2 - 6t
            dH10 = 3t^2 - 4t + 1
            dH01 = -6t^2 + 6t
            dH11 = 3t^2 - 2t
            tx = dH00*p0[1] + dH10*m0s[1] + dH01*p1[1] + dH11*m1s[1]
            ty = dH00*p0[2] + dH10*m0s[2] + dH01*p1[2] + dH11*m1s[2]
            tz = dH00*p0[3] + dH10*m0s[3] + dH01*p1[3] + dH11*m1s[3]
            tunit = normalize3((tx,ty,tz))
            push!(P, (px,py,pz))
            push!(T, tunit)
            push!(seg_of, i)
        end
    end
    return P, T, seg_of
end

# --------- Parallel-transport frames along a space curve ----------
"""
    build_frames(P, T)

Parallel-transport frames (N, B) minimizing twist: returns two arrays `N`, `B`
(each element is a unit vector), with `length(P) == length(N) == length(B)`.
"""
function build_frames(P::Vector{NTuple{3,Float64}},
                    T::Vector{NTuple{3,Float64}})
    N = NTuple{3,Float64}[]
    B = NTuple{3,Float64}[]

    # Initial frame: pick an "up" not parallel to T[1]
    up_candidates = ((0.0,0.0,1.0), (0.0,1.0,0.0), (1.0,0.0,0.0))
    up = abs(dot(T[1], up_candidates[1])) < 0.9 ? up_candidates[1] :
        abs(dot(T[1], up_candidates[2])) < 0.9 ? up_candidates[2] : up_candidates[3]
    n0 = normalize3(cross(up, T[1]))
    b0 = cross(T[1], n0)
    push!(N, n0); push!(B, b0)

    for i in 2:length(P)
        Ti_1, Ti = T[i-1], T[i]
        k = cross(Ti_1, Ti)
        s = norm3(k)
        if s < 1e-12
            # Nearly parallel: re-orthonormalize previous normal against current tangent
            nprev = N[i-1]
            nproj = nprev .- (dot(nprev, Ti)) .* Ti
            ni = normalize3((nproj[1], nproj[2], nproj[3]))
            bi = cross(Ti, ni)
            push!(N, ni); push!(B, bi)
        else
            k̂ = (k[1]/s, k[2]/s, k[3]/s)
            θ  = atan(s, dot(Ti_1, Ti)) # atan2(||cross||, dot)
            ni = rotate_around_axis(N[i-1], k̂, θ)
            bi = cross(Ti, ni)
            push!(N, ni); push!(B, bi)
        end
    end
    return N, B
end

# --------- Ring vertex generation ----------
function build_rings(P, N, B; radius::Real, ntheta::Int)
    verts = NTuple{3,Float64}[]
    vert_cell = Int[]
    ring_offset = Int[1]  # 1-based start of each ring in 'verts'
    for i in eachindex(P)
        for j in 0:(ntheta-1)
            θ = 2π * j/ntheta
            nx, ny, nz = N[i]
            bx, by, bz = B[i]
            cx = radius*(cos(θ)*nx + sin(θ)*bx)
            cy = radius*(cos(θ)*ny + sin(θ)*by)
            cz = radius*(cos(θ)*nz + sin(θ)*bz)
            px = P[i][1] + cx
            py = P[i][2] + cy
            pz = P[i][3] + cz
            push!(verts, (px,py,pz))
            push!(vert_cell, i)  # ring index as "cell" owner for vertices
        end
        push!(ring_offset, length(verts)+1)
    end
    return verts, vert_cell, ring_offset
end

# --------- Tube sidewall triangulation (between consecutive rings) ----------
function stitch_sidewalls(nrings::Int, ntheta::Int, ring_offset::Vector{Int})
    tris = Vector{NTuple{3,Int}}()
    face_cell = Int[]  # segment index 1..(nrings-1)
    for i in 1:(nrings-1)
        baseA = ring_offset[i]
        baseB = ring_offset[i+1]
        for j in 0:(ntheta-1)
            a0 = baseA + j
            a1 = baseA + ((j+1) % ntheta)
            b0 = baseB + j
            b1 = baseB + ((j+1) % ntheta)
            # Two triangles per quad (consistent winding: outward normal)
            push!(tris, (a0, a1, b1))
            push!(tris, (a0, b1, b0))
            push!(face_cell, i)
            push!(face_cell, i)
        end
    end
    return tris, face_cell
end

# --------- Caps: flat or rounded (hemisphere) ----------
function cap_flat!(side_verts::Vector{NTuple{3,Float64}},
                side_vert_cell::Vector{Int},
                ring_offset::Vector{Int},
                on_start::Bool, ntheta::Int)
    # center vertex index
    ring = on_start ? 1 : (length(ring_offset)-1)
    first = ring_offset[ring]
    last  = first + ntheta - 1
    # center = side_verts[first:first] |> first
    # We’ll compute exact center from average of ring vertices (more robust)
    cx = 0.0; cy = 0.0; cz = 0.0
    for k in first:last
        cx += side_verts[k][1]; cy += side_verts[k][2]; cz += side_verts[k][3]
    end
    cx /= ntheta; cy /= ntheta; cz /= ntheta
    push!(side_verts, (cx,cy,cz))
    push!(side_vert_cell, ring)

    center_idx = length(side_verts)
    cap_tris   = Vector{NTuple{3,Int}}()
    cap_cells  = Int[]

    if on_start
        # fan: (j+1, j, center)
        for j in 0:(ntheta-1)
            v0 = first + ((j+1) % ntheta)
            v1 = first + j
            push!(cap_tris, (v0, v1, center_idx))
            push!(cap_cells, ring)
        end
    else
        # fan: (j, j+1, center) to keep outward normals
        for j in 0:(ntheta-1)
            v0 = first + j
            v1 = first + ((j+1) % ntheta)
            push!(cap_tris, (v0, v1, center_idx))
            push!(cap_cells, ring-1)  # associate with last segment
        end
    end
    return cap_tris, cap_cells
end

function cap_hemisphere!(P::Vector{NTuple{3,Float64}},
                        T::Vector{NTuple{3,Float64}},
                        N::Vector{NTuple{3,Float64}},
                        B::Vector{NTuple{3,Float64}},
                        side_verts::Vector{NTuple{3,Float64}},
                        side_vert_cell::Vector{Int},
                        ring_offset::Vector{Int};
                        radius::Real, ntheta::Int, stacks::Int=4, on_start::Bool)
    # Choose anchor ring and local frame
    idx = on_start ? 1 : length(P)
    # Hemisphere center along -T (start) or +T (end)
    dir = on_start ? (-T[idx][1], -T[idx][2], -T[idx][3]) : T[idx]
    c0  = (P[idx][1] + radius*dir[1],
        P[idx][2] + radius*dir[2],
        P[idx][3] + radius*dir[3])

    # Build additional rings (excluding the boundary ring we already have)
    new_ring_offsets = Int[]
    for s in 1:stacks
        φ = (π/2) * s/stacks
        r_ring = radius * cos(φ)      # shrinking radius
        # shift along ±T so the hemisphere meets the side ring
        offset = radius * (on_start ? -sin(φ) : sin(φ))
        center = (P[idx][1] + offset*T[idx][1],
                P[idx][2] + offset*T[idx][2],
                P[idx][3] + offset*T[idx][3])
        push!(new_ring_offsets, length(side_verts)+1)
        for j in 0:(ntheta-1)
            θ = 2π * j/ntheta
            nx, ny, nz = N[idx]
            bx, by, bz = B[idx]
            px = center[1] + r_ring*(cos(θ)*nx + sin(θ)*bx)
            py = center[2] + r_ring*(cos(θ)*ny + sin(θ)*by)
            pz = center[3] + r_ring*(cos(θ)*nz + sin(θ)*bz)
            push!(side_verts, (px,py,pz))
            push!(side_vert_cell, on_start ? 1 : (length(P)-1))
        end
    end

    # Triangulate stacks as sidewalls (between successive rings),
    # plus the final tip as a fan
    cap_tris  = Vector{NTuple{3,Int}}()
    cap_cells = Int[]

    # From boundary ring to first inner ring
    baseBoundary = ring_offset[on_start ? 1 : (length(ring_offset)-1)]
    prev = baseBoundary
    # iterate rings inwards
    ringsA = vcat(prev, new_ring_offsets[1:end-1])
    ringsB = new_ring_offsets

    for (rA, rB) in zip(ringsA, ringsB)
        for j in 0:(ntheta-1)
            a0 = rA + j
            a1 = rA + ((j+1) % ntheta)
            b0 = rB + j
            b1 = rB + ((j+1) % ntheta)
            if on_start
                push!(cap_tris, (a1, a0, b1))
                push!(cap_tris, (a0, b0, b1))
            else
                push!(cap_tris, (a0, a1, b1))
                push!(cap_tris, (a0, b1, b0))
            end
            push!(cap_cells, on_start ? 1 : (length(P)-1))
            push!(cap_cells, on_start ? 1 : (length(P)-1))
        end
    end

    # Tip vertex (cap center)
    tip = (c0[1], c0[2], c0[3])
    push!(side_verts, tip)
    push!(side_vert_cell, on_start ? 1 : (length(P)-1))
    tip_idx = length(side_verts)
    last_ring = new_ring_offsets[end]

    if on_start
        for j in 0:(ntheta-1)
            v0 = last_ring + ((j+1) % ntheta)
            v1 = last_ring + j
            push!(cap_tris, (v0, v1, tip_idx))
            push!(cap_cells, 1)
        end
    else
        for j in 0:(ntheta-1)
            v0 = last_ring + j
            v1 = last_ring + ((j+1) % ntheta)
            push!(cap_tris, (v0, v1, tip_idx))
            push!(cap_cells, length(P)-1)
        end
    end

    return cap_tris, cap_cells
end

"""
    triangulate_cylindrical_mesh(points;
        radius, ntheta=24, samples_per_seg=8,
        cap=:flat, cap_stacks=4)

Triangulate a *tube* (cylindrical surface) around a polyline given by `points`
(2×M or 3×M). The polyline is smoothed with a Catmull–Rom spline.

# Keyword arguments
- `radius::Real`          : Tube radius (required).
- `ntheta::Int=24`        : Vertices per ring (circumferential resolution).
- `samples_per_seg::Int=8`: Spline samples per original segment.
- `cap::Symbol=:flat`     : `:flat` (disks) or `:round` (hemispherical).
- `cap_stacks::Int=4`     : Number of stacks for rounded caps.

# Returns
NamedTuple with:
- `xyz::Matrix{Float64}`      (3×N)
- `tri::Matrix{Int}`          (3×NT)
- `face_cell::Vector{Int}`    (NT)
- `vert_cell::Vector{Int}`    (N)
- `ring_offset::Vector{Int}`  (nrings+1)

The `*_cell` arrays assign each face/vertex to a “cell” (segment index) along
the tube, so you can color/select by cell.
"""
function triangulate_cylindrical_mesh(points::AbstractMatrix{<:Real};
        radius::Real, ntheta::Int=24, samples_per_seg::Int=8,
        cap::Symbol=:flat, cap_stacks::Int=4)

    @assert size(points,2) ≥ 2 "Need at least two points"
    @assert radius > 0 "radius must be positive"
    @assert ntheta ≥ 3 "ntheta must be ≥ 3"
    @assert samples_per_seg ≥ 1 "samples_per_seg must be ≥ 1"
    @assert cap in (:flat, :round) "cap must be :flat or :round"

    # 1) Smooth centerline and tangents
    P, T, seg_of = sample_centerline(points; samples_per_seg)

    # 2) Frames (parallel transport)
    N, B = build_frames(P, T)

    # 3) Rings + sidewall vertices
    side_verts, vert_cell, ring_offset = build_rings(P, N, B; radius, ntheta)
    nrings = length(P)

    # 4) Sidewall triangles
    side_tris, face_cell = stitch_sidewalls(nrings, ntheta, ring_offset)

    # 5) Caps
    cap_tris = NTuple{3,Int}[]
    cap_cells = Int[]
    if cap == :flat
        t1, c1 = cap_flat!(side_verts, vert_cell, ring_offset, true, ntheta)
        t2, c2 = cap_flat!(side_verts, vert_cell, ring_offset, false, ntheta)
        append!(cap_tris, t1); append!(cap_cells, c1)
        append!(cap_tris, t2); append!(cap_cells, c2)
    else
        t1, c1 = cap_hemisphere!(P, T, N, B, side_verts, vert_cell, ring_offset; radius, ntheta, stacks=cap_stacks, on_start=true)
        t2, c2 = cap_hemisphere!(P, T, N, B, side_verts, vert_cell, ring_offset; radius, ntheta, stacks=cap_stacks, on_start=false)
        append!(cap_tris, t1); append!(cap_cells, c1)
        append!(cap_tris, t2); append!(cap_cells, c2)
    end

    # 6) Assemble outputs
    all_tris    = vcat(side_tris, cap_tris)
    all_cells   = vcat(face_cell, cap_cells)

    # xyz as 3×N Matrix
    Nvert = length(side_verts)
    xyz = Matrix{Float64}(undef, Nvert, 3)
    for (i,v) in enumerate(side_verts)
        xyz[i, 1] = v[1]
        xyz[i, 2] = v[2]
        xyz[i, 3] = v[3]
    end
    # tri as 3×NT Matrix
    NTri = length(all_tris)
    tri = Matrix{Int}(undef, NTri, 3)
    for (i,(a,b,c)) in enumerate(all_tris)
        tri[i, 1] = a;
        tri[i, 2] = b;
        tri[i, 3] = c
    end

    return (xyz = xyz,
            tri = tri,
            face_cell = all_cells,
            vert_cell = vert_cell,
            ring_offset = ring_offset)
end