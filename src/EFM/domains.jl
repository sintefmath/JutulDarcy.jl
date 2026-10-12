function normal_conductivity(values, cell, normal)
    if values isa AbstractVector
        return values[cell]
    end
    K = Jutul.expand_perm(view(values, :, cell), 3)
    return dot(normal, K*normal)
end

function series_conductivity(area, dl, dr, kl, kr)
    if kl == 0 || kr == 0
        return 0.0
    end
    return area/(dl/kl + dr/kr)
end

function fracture_conductivity(domain, values)
    mesh = physical_representation(domain)::EmbeddedFractureMesh
    N = tpfv_geometry(mesh).neighbors
    aperture = domain[:aperture]
    T = zeros(Jutul.float_type(mesh), size(N, 2))
    for i in eachindex(T)
        l, r = N[:, i]
        dl, dr = mesh.geometry.connection_distances[:, i]
        T[i] = series_conductivity(tpfv_geometry(mesh).areas[i], dl, dr, values[l]*aperture[l], values[r]*aperture[r])
    end
    for intersection in mesh.geometry.intersections
        cells = intersection.cells
        half = intersection.length.*values[cells].*aperture[cells]./intersection.distances
        denominator = sum(half)
        for face in intersection.faces
            l, r = N[:, face]
            il, ir = findfirst(==(l), cells), findfirst(==(r), cells)
            T[face] = 0.0
            if denominator > 0
                T[face] = half[il]*half[ir]/denominator
            end
        end
    end
    return T
end

function matrix_fracture_conductivity(domain, matrix_values, fracture_values)
    mesh = physical_representation(domain)::EmbeddedFractureMesh
    aperture = domain[:aperture]
    T = zeros(Jutul.float_type(mesh), length(mesh.geometry.matrix_cells))
    for i in eachindex(T)
        mc, fc = mesh.geometry.matrix_cells[i], mesh.geometry.fracture_cells[i]
        normal = mesh.geometry.frames[mesh.geometry.patches[fc].fracture].normal
        km = normal_conductivity(matrix_values, mc, normal)
        kf = fracture_values[fc]
        T[i] = series_conductivity(mesh.geometry.matrix_areas[i], mesh.geometry.matrix_distances[i], aperture[fc]/2, km, kf)
    end
    return T
end

function fracture_property(values, default, ids, name)
    T = eltype(default)
    if isnothing(values)
        out = default
    elseif values isa Number
        out = fill(T(values), length(ids))
    else
        # Vector overrides are per fracture object, not per clipped patch.
        length(values) >= maximum(ids) || throw(ArgumentError("$name must contain one value per fracture"))
        out = T.(values[ids])
    end
    all(isfinite, out) || throw(ArgumentError("$name must be finite"))
    return out
end

"""
    embedded_fracture_domain(matrix, fractures; aperture=nothing, permeability=nothing, porosity=nothing, kwarg...)
    embedded_fracture_domain(fractures, matrix; kwarg...)

Create a fracture `DataDomain` with clipped EFM cells, fracture–fracture
transmissibilities, and DFM-compatible `FractureMatrixConnection` data. The
matrix `DataDomain` and its mesh remain unchanged. Physical properties default
to the fracture objects; scalar or per-fracture keyword overrides are accepted.
Other keywords are passed to `reservoir_domain`. Geometry options are
`tolerance`, `intersections`, and `float_type` (defaults to the matrix mesh's
floating-point type). Permeability uses SI units and can be a full
symmetric tensor in the matrix. Fracture permeability is isotropic in-plane.

The standard EDFM approximation is intended for conductive fractures; interior
matrix connections remain open, so impermeable fracture barriers require a different
discretization such as DFM or pEDFM.
"""
function embedded_fracture_domain(matrix::DataDomain, fractures;
        tolerance = nothing,
        intersections = true,
        float_type = Jutul.float_type(physical_representation(matrix)),
        kwarg...
    )
    mesh = embedded_fracture_mesh(physical_representation(matrix), fractures; tolerance = tolerance, intersections = intersections, float_type = float_type)
    return embedded_fracture_domain(mesh, matrix; kwarg...)
end

embedded_fracture_domain(fractures::Union{AbstractEmbeddedFracture, AbstractVector, Tuple}, matrix::DataDomain; kwarg...) = embedded_fracture_domain(matrix, fractures; kwarg...)

function embedded_fracture_domain(mesh::EmbeddedFractureMesh, matrix::DataDomain;
        aperture = nothing,
        hydraulic_aperture = nothing,
        permeability = nothing,
        porosity = nothing,
        kwarg...
    )
    physical_representation(matrix) === mesh.matrix_mesh || throw(ArgumentError("fracture mesh must be used with its original matrix mesh"))
    T = Jutul.float_type(mesh)
    ids = [p.fracture for p in mesh.geometry.patches]
    properties = fracture_properties(mesh.fractures, ids, T; aperture = aperture,
        hydraulic_aperture = hydraulic_aperture, permeability = permeability, porosity = porosity)
    a, k, phi = properties.aperture, properties.permeability, properties.porosity
    domain = reservoir_domain(mesh; permeability = k, porosity = phi, kwarg...)
    domain[:aperture, Cells()] = a
    domain[:fracture_id, Cells()] = ids
    domain[:matrix_face_coverage, nothing] = matrix_face_coverage(mesh)
    domain[:volumes, Cells()] = tpfv_geometry(mesh).volumes.*a
    N = tpfv_geometry(mesh).neighbors
    domain[:areas, Faces()] = tpfv_geometry(mesh).areas.*vec(sum(a[N]; dims = 1))/2
    domain[:cell_normals, Cells()] = reduce(hcat, [mesh.geometry.frames[i].normal for i in ids])
    domain[:transmissibilities, Faces()] = fracture_conductivity(domain, k)
    domain[:rock_thermal_conductivities, Faces()] = fracture_conductivity(domain, vec(domain[:rock_thermal_conductivity]).*(1 .- phi))
    fluid = domain[:fluid_thermal_conductivity]
    if fluid isa AbstractVector
        domain[:fluid_thermal_conductivities, Faces()] = fracture_conductivity(domain, fluid.*phi)
    else
        domain[:fluid_thermal_conductivities, Faces()] = reduce(vcat, [fracture_conductivity(domain, fluid[ph, :].*phi)' for ph in axes(fluid, 1)])
    end
    fmc = FractureMatrixConnection()
    mc, fc = mesh.geometry.matrix_cells, mesh.geometry.fracture_cells
    domain.entities[fmc] = length(mc)
    domain[:matrix_cells, fmc] = mc
    domain[:connection_cells, fmc] = fc
    domain[:matrix_cell_centroids, fmc] = matrix[:cell_centroids][:, mc]
    domain[:matrix_distance, fmc] = mesh.geometry.matrix_distances
    domain[:matrix_connection_area, fmc] = mesh.geometry.matrix_areas
    domain[:fracture_matrix_transmissibility, fmc] = matrix_fracture_conductivity(domain, matrix[:permeability], k)
    domain[:fracture_matrix_gravity_difference, fmc] = JutulDarcy.gravity_constant*(matrix[:cell_centroids][3, mc] .- domain[:cell_centroids][3, fc])
    # Match the DFM thermal cross term, with rock and pore-fluid contributions.
    mr = vec(matrix[:rock_thermal_conductivity]).*(1 .- vec(matrix[:porosity]))
    fr = vec(domain[:rock_thermal_conductivity]).*(1 .- phi)
    mf = matrix[:fluid_thermal_conductivity]
    ff = domain[:fluid_thermal_conductivity]
    if (mf isa AbstractMatrix && size(mf, 1) != 1) || (ff isa AbstractMatrix && size(ff, 1) != 1)
        throw(ArgumentError("EFM cross-term thermal conductivity currently requires one fluid conductivity per cell"))
    end
    thermal = matrix_fracture_conductivity(domain, mr, fr)
    thermal += matrix_fracture_conductivity(domain, vec(mf).*vec(matrix[:porosity]), vec(ff).*phi)
    domain[:fracture_matrix_thermal_conductivity, fmc] = thermal
    return domain
end

function matrix_face_coverage(mesh::EmbeddedFractureMesh)
    matrix = UnstructuredMesh(mesh.matrix_mesh)
    covered = zeros(number_of_faces(matrix))
    face_areas = tpfv_geometry(matrix).areas
    for patch in mesh.geometry.patches
        length(patch.matrix_cells) > 1 || continue
        frame = mesh.geometry.frames[patch.fracture]
        for face in coincident_faces(matrix, first(patch.matrix_cells), frame, mesh.geometry.tolerance)
            all(c -> c in patch.matrix_cells, matrix.faces.neighbors[face]) || continue
            points = [project(frame, matrix.node_points[n]) for n in matrix.faces.faces_to_nodes[face]]
            area = 0.0
            for boundary in convex_parts(points, mesh.geometry.tolerance), part in patch.parts
                clipped = clip_polygon(part, boundary, mesh.geometry.tolerance)
                if length(clipped) >= 3
                    area += polygon_area_centroid(clipped)[1]
                end
            end
            covered[face] += area/face_areas[face]
        end
    end
    return clamp.(covered, 0.0, 1.0)
end

function embedded_matrix_domain(matrix, fractures; subtract_fracture_volume = false)
    coverage = fractures[:matrix_face_coverage]
    if !subtract_fracture_volume && all(==(0), coverage)
        return matrix
    end
    adjusted = DataDomain(physical_representation(matrix))
    for (key, (value, entity)) in pairs(matrix)
        adjusted[key, entity] = deepcopy(value)
    end
    if subtract_fracture_volume
        subtract_volumes!(adjusted[:volumes], fractures)
    end
    if any(>(0), coverage)
        # A face-coincident fracture replaces the covered part of the direct
        # matrix connection. Keeping that path would count normal flow twice.
        model = reservoir_model(setup_reservoir_model(matrix, SinglePhaseSystem(AqueousPhase())))
        variables = ((:transmissibilities, JutulDarcy.Transmissibilities(), :Transmissibilities),
            (:rock_thermal_conductivities, JutulDarcy.RockThermalConductivities(), :RockThermalConductivities),
            (:fluid_thermal_conductivities, JutulDarcy.FluidThermalConductivities(), :FluidThermalConductivities))
        for (key, variable, name) in variables
            values = Jutul.default_parameter_values(matrix, model, variable, name)
            if values isa AbstractVector
                values .*= 1 .- coverage
            else
                values .*= (1 .- coverage)'
            end
            adjusted[key, Faces()] = values
        end
    end
    return adjusted
end

"""
    combined_fracture_domain(matrix, fractures; subtract_fracture_volume=false)

Combine an EFM fracture domain with the matrix into one TPFA connection graph.
The original matrix cell numbers and faces come first; fracture cells and NNCs
are appended. Matrix boundary geometry is retained. Explicit transmissibilities
are used for all added connections. The original matrix mesh is unchanged.
The integer cell field `:is_fracture` is `0` for matrix cells and `1` for fracture
cells.
"""
function combined_fracture_domain(matrix::DataDomain, fractures::DataDomain; subtract_fracture_volume = false)
    fm = physical_representation(fractures)::EmbeddedFractureMesh
    physical_representation(matrix) === fm.matrix_mesh || throw(ArgumentError("EFM domains must share the original matrix mesh"))
    matrix = embedded_matrix_domain(matrix, fractures)
    nm = number_of_cells(matrix)
    fmc = FractureMatrixConnection()
    fc = fractures[:connection_cells, fmc]
    volumes = vcat(matrix[:volumes], fractures[:volumes])
    if subtract_fracture_volume
        subtract_volumes!(view(volumes, 1:nm), fractures)
    end
    domain = unified_fracture_domain(matrix, fractures, fm.geometry.matrix_areas,
        fractures[:cell_normals][:, fc]; volumes = volumes)
    # Use the matrix's standard parameter initialization so overrides,
    # multipliers, and existing NNCs are incorporated exactly once.
    mm = reservoir_model(setup_reservoir_model(matrix, SinglePhaseSystem(AqueousPhase())))
    Tm = Jutul.default_parameter_values(matrix, mm, JutulDarcy.Transmissibilities(), :Transmissibilities)
    domain[:transmissibilities, Faces()] = vcat(Tm, fractures[:transmissibilities], fractures[:fracture_matrix_transmissibility, fmc])
    Tr = Jutul.default_parameter_values(matrix, mm, JutulDarcy.RockThermalConductivities(), :RockThermalConductivities)
    Tf = Jutul.default_parameter_values(matrix, mm, JutulDarcy.FluidThermalConductivities(), :FluidThermalConductivities)
    mr = vec(matrix[:rock_thermal_conductivity]).*(1 .- vec(matrix[:porosity]))
    fr = vec(fractures[:rock_thermal_conductivity]).*(1 .- vec(fractures[:porosity]))
    cross_rock = matrix_fracture_conductivity(fractures, mr, fr)
    cross_fluid = fractures[:fracture_matrix_thermal_conductivity, fmc]-cross_rock
    domain[:rock_thermal_conductivities, Faces()] = vcat(Tr, fractures[:rock_thermal_conductivities], cross_rock)
    domain[:fluid_thermal_conductivities, Faces()] = vcat(vec(Tf), vec(fractures[:fluid_thermal_conductivities]), cross_fluid)
    return domain
end

function subtract_volumes!(volumes, fractures)
    mesh = physical_representation(fractures)::EmbeddedFractureMesh
    removed = zeros(length(volumes))
    for (cell, patch) in enumerate(mesh.geometry.patches)
        for mc in patch.matrix_cells
            removed[mc] += fractures[:volumes][cell]/length(patch.matrix_cells)
        end
    end
    all(volumes .> removed) || throw(ArgumentError("fracture volume exceeds an intersected matrix cell's volume"))
    volumes .-= removed
    return volumes
end
