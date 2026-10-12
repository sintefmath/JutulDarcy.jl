using Jutul, JutulDarcy, Test, LinearAlgebra, SparseArrays
import JutulDarcy.EFM

function efm_connection_matrix(domain, T, cell_map = collect(1:number_of_cells(domain)))
    N = domain[:neighbors]
    return sparse(cell_map[N[1, :]], cell_map[N[2, :]], T, length(cell_map), length(cell_map)) +
           sparse(cell_map[N[2, :]], cell_map[N[1, :]], T, length(cell_map), length(cell_map))
end

function efm_dfm_geometry(; cross = false, three = false, n = (4, 4, 4), L = (4.0, 4.0, 4.0), aperture = 1e-4)
    mesh = reservoir_mesh(n, L)
    matrix = reservoir_domain(mesh; permeability = 1e-2si_unit(:darcy), porosity = 0.05)
    objects = [RectangleFracture(L./2, [1.0, 0.0, 0.0], (L[2], L[3]); aperture = aperture)]
    normals = [[1.0, 0.0, 0.0]]
    if cross
        push!(objects, RectangleFracture(L./2, [0.0, 1.0, 0.0], (L[1], L[3]); aperture = aperture))
        push!(normals, [0.0, 1.0, 0.0])
    end
    if three
        push!(objects, RectangleFracture(L./2, [0.0, 0.0, 1.0], (L[1], L[2]); aperture = aperture))
        push!(normals, [0.0, 0.0, 1.0])
    end
    cuts = [Jutul.PlaneCut(collect(L./2), normal) for normal in normals]
    cut, info = Jutul.cut_mesh(mesh, cuts; extra_out = true)
    faces = findall(==(0), info[:face_index])
    # Aligned cuts reuse existing faces rather than creating new ones.
    if isempty(faces)
        geo = tpfv_geometry(cut)
        for (i, normal) in enumerate(normals)
            axis = argmax(abs.(normal))
            append!(faces, findall(f -> abs(geo.face_centroids[axis, f]-L[axis]/2) < 1e-10 && abs(geo.normals[axis, f]) > 0.99, 1:number_of_faces(cut)))
        end
        unique!(faces)
    end
    fmesh = Jutul.EmbeddedMesh(cut, faces; intersection_strategy = :star_delta)
    dfm_matrix = reservoir_domain(cut; permeability = 1e-2si_unit(:darcy), porosity = 0.05)
    dfm = fracture_domain(fmesh, dfm_matrix; aperture = aperture, porosity = 0.5)
    efm = embedded_fracture_domain(matrix, objects)
    return (; matrix, efm, dfm_matrix, dfm, info)
end

function efm_validation_case(matrix, fractures, system;
        representation = :coupled,
        thermal = false,
        nstep = 8,
        total_time = 2000.0,
        rate = 1e-5,
        dfm = false,
        well_cells = nothing,
        subtract_fracture_volume = false
    )
    mesh = physical_representation(matrix)
    if isnothing(well_cells)
        dims = Jutul.grid_dims_ijk(mesh)
        inj_cells = findall(c -> cell_ijk(mesh, c)[1:2] == (1, 1), 1:number_of_cells(mesh))
        prod_cells = findall(c -> cell_ijk(mesh, c)[1:2] == dims[1:2], 1:number_of_cells(mesh))
    else
        inj_cells, prod_cells = well_cells
    end
    inj = setup_well(matrix, inj_cells; name = :Injector, simple_well = true)
    prod = setup_well(matrix, prod_cells; name = :Producer, simple_well = true)
    if dfm
        model = setup_fractured_reservoir_model(matrix, fractures, system; wells = [inj, prod], thermal = thermal)
    else
        model = setup_embedded_fracture_model(matrix, fractures, system;
            wells = [inj, prod], thermal = thermal, representation = representation,
            subtract_fracture_volume = subtract_fracture_volume)
    end
    mix = [1.0]
    initial = (Pressure = 50si_unit(:bar),)
    if system isa ImmiscibleSystem && number_of_phases(system) == 2
        mix = [1.0, 0.0]
        initial = (initial..., Saturations = [0.0, 1.0])
    end
    injection = (density = 1000.0,)
    if system isa ImmiscibleSystem
        injection = (density = system.rho_ref[1],)
    end
    if thermal
        initial = (initial..., Temperature = 363.15)
        injection = (injection..., temperature = 283.15)
    end
    controls = Dict(
        :Injector => InjectorControl(TotalRateTarget(rate), mix; injection...),
        :Producer => ProducerControl(BottomHolePressureTarget(10si_unit(:bar))))
    state = setup_reservoir_state(model; initial...)
    forces = setup_reservoir_forces(model; control = controls)
    return JutulCase(model, fill(total_time/nstep, nstep), forces; state0 = state)
end

function efm_dfm_comparison(n; thermal = false, pv_fraction = 0.5)
    L = (10.0, 10.0, 3.0)
    d = efm_dfm_geometry(; cross = true, three = true, n = n, L = L)
    system = ImmiscibleSystem((AqueousPhase(), LiquidPhase()); reference_densities = (800.0, 1000.0))
    if thermal
        system = :geothermal
    end
    mesh = physical_representation(d.matrix)
    ijk = [cell_ijk(mesh, c) for c in 1:number_of_cells(mesh)]
    mapping = d.info[:cell_index]
    inj = findall(c -> ijk[mapping[c]][1:2] == (1, 1), eachindex(mapping))
    prod = findall(c -> ijk[mapping[c]][1:2] == n[1:2], eachindex(mapping))
    total_time = 0.1si_unit(:year)
    rate = pv_fraction*(sum(pore_volume(d.matrix))+sum(pore_volume(d.efm)))/total_time
    efm = efm_validation_case(d.matrix, d.efm, system; thermal = thermal,
        total_time = total_time, nstep = 20, rate = rate, subtract_fracture_volume = true)
    dfm = efm_validation_case(d.dfm_matrix, d.dfm, system; thermal = thermal,
        total_time = total_time, nstep = 20, rate = rate, dfm = true, well_cells = (inj, prod))
    solver = GenericKrylov(preconditioner = ILUZeroPreconditioner(), relative_tolerance = 1e-10)
    simarg = (info_level = -1, linear_solver = solver, tol_cnv = 1e-7,
        tol_cnve = 1e-7, initial_dt = 5.0, max_nonlinear_iterations = 20)
    re = simulate_reservoir(efm; simarg...)
    rd = simulate_reservoir(dfm; simarg..., linear_solver = GenericKrylov(preconditioner = ILUZeroPreconditioner(), relative_tolerance = 1e-10))
    # DFM divides matrix cells; compare its volume-weighted averages
    # with EFM's original matrix cells at the same reporting times.
    volumes = d.dfm_matrix[:volumes]
    weights = zeros(number_of_cells(d.matrix))
    for (i, c) in enumerate(mapping)
        weights[c] += volumes[i]
    end
    function average_dfm(values)
        out = zeros(length(weights))
        for (i, c) in enumerate(mapping)
            out[c] += values[i]*volumes[i]
        end
        return out./weights
    end
    pressure_error, transport_error = 0.0, 0.0
    for (e, r) in zip(re.result.states, rd.result.states)
        pe = e[:Reservoir][:Pressure]
        pd = average_dfm(r[:Reservoir][:Pressure])
        pressure_error = max(pressure_error, sqrt(sum(((pe-pd)./pd).^2)/length(pd)))
        if thermal
            te = e[:Reservoir][:Temperature]
            td = average_dfm(r[:Reservoir][:Temperature])
        else
            te = e[:Reservoir][:Saturations][1, :]
            td = average_dfm(r[:Reservoir][:Saturations][1, :])
        end
        transport_error = max(transport_error, sqrt(sum((te-td).^2)/length(td)))
    end
    return (; pressure_error, transport_error, efm = re, dfm = rd)
end

@testset "EFM" begin
    @testset "Fracture objects and exact clipping" begin
        mesh = reservoir_mesh((3, 3, 3), (3.0, 3.0, 3.0))
        matrix = reservoir_domain(mesh)
        before = deepcopy(matrix[:volumes])
        coordinates = copy(mesh.node_points)
        f = RectangleFracture([1.5, 1.5, 1.5], [1, 0, 0], (3, 3))
        domain = embedded_fracture_domain(matrix, f)
        @test number_of_cells(domain) == 9
        @test number_of_faces(domain) == 12
        @test sum(domain[:volumes]) ≈ 9f.properties.aperture
        @test domain[:matrix_distance] ≈ fill(0.25, 9)
        @test physical_representation(domain).matrix_mesh === mesh
        @test matrix[:volumes] == before
        @test mesh.node_points == coordinates
        @test all(>(0), domain[:fracture_matrix_transmissibility])
        @test all(>(0), domain[:transmissibilities])
        @test embedded_fracture_domain(f, matrix)[:transmissibilities] == domain[:transmissibilities]
        @test all(isfinite, setup_parameters(setup_embedded_fracture_model(matrix, domain, SinglePhaseSystem()))[:Fractures][:Transmissibilities])

        ellipse = EllipticalFracture([1.5, 1.5, 1.5], [1, 0, 0], (1.0, 0.75); segments = 128)
        ed = embedded_fracture_domain(matrix, ellipse)
        @test sum(ed[:volumes])/ellipse.properties.aperture ≈ pi*0.75 rtol = 5e-4
        @test EllipsoidFracture === EllipticalFracture
        polygon = PolygonFracture([[1.5, 0, 0], [1.5, 2, 0], [1.5, 2, 1], [1.5, 1, 1], [1.5, 1, 2], [1.5, 0, 2]])
        pd = embedded_fracture_domain(matrix, polygon)
        @test sum(pd[:volumes])/polygon.properties.aperture ≈ 3
        # Clipping a concave surface can leave disconnected pieces inside one
        # cell. These must remain distinct fracture cells.
        one = reservoir_domain(reservoir_mesh((1, 1, 1), (1.0, 1.0, 1.0)))
        disconnected = PolygonFracture([[0.5, -1, 0.1], [0.5, 0.8, 0.1], [0.5, 0.8, 0.3],
            [0.5, -0.5, 0.3], [0.5, -0.5, 0.7], [0.5, 0.8, 0.7], [0.5, 0.8, 0.9], [0.5, -1, 0.9]])
        disconnected_domain = embedded_fracture_domain(one, disconnected)
        @test number_of_cells(disconnected_domain) == 2
        @test number_of_faces(disconnected_domain) == 0
        @test sum(disconnected_domain[:volumes])/disconnected.properties.aperture ≈ 0.32
        rotated = RectangleFracture([1.5, 1.5, 1.5], [1, 1, 1], (0.5, 0.8); angle = 0.3)
        rd = embedded_fracture_domain(matrix, rotated)
        @test sum(rd[:volumes])/rotated.properties.aperture ≈ 0.4
        # Non-convex face fans can retrace a collinear boundary segment. The
        # retraced portion has zero signed area and must not obstruct clipping.
        retraced = EFM.Point2.([(0, 0), (3, 0), (2, 0), (4, 0), (4, 1), (0, 1)])
        @test sum(p -> EFM.polygon_area_centroid(p)[1], EFM.convex_parts(retraced, 1e-9)) ≈ 4
        crossing = EFM.Point2.([(0, 0), (2, 0), (2, 2), (0, 2), (0, 0),
            (1, 1), (3, 1), (3, 3), (1, 3), (1, 1), (0, 0)])
        @test sum(p -> EFM.polygon_area_centroid(p)[1], EFM.convex_parts(crossing, 1e-9; resolve_winding = true)) ≈ 7
        sliver = RectangleFracture([0.5, 0.5, 0.5], [1, 0, 0], (1e-10, 1.0))
        @test_throws ArgumentError embedded_fracture_domain(one, sliver; tolerance = 1e-8)
        # Exact distance integration includes tetrahedra split two-and-two.
        p = (EFM.Point3(0, 0, 0), EFM.Point3(1, 0, 0), EFM.Point3(0, 1, 0), EFM.Point3(0, 0, 1))
        @test EFM.tetra_average_distance(p, EFM.SVector(1.0, 1.0, -1.0, -1.0)) ≈ 3/8
        @test EFM.tetra_average_distance(p, EFM.SVector(1.0, -1.0, -1.0, -1.0)) ≈ 9/16
        @test_throws ArgumentError PolygonFracture([[0, 0, 0], [1, 0, 0], [1, 1, 1], [0, 1, 0]])
        @test_throws ArgumentError PolygonFracture([[0, 0, 0], [2, 0, 0], [0, 1, 0], [2, 2, 0]])
        @test_throws ArgumentError RectangleFracture([0, 0, 0], [0, 0, 0], (1, 1))
        @test_throws ArgumentError EllipticalFracture([0, 0, 0], [0, 0, 1], (-1, 1))
        @test_throws ArgumentError embedded_fracture_domain(matrix, RectangleFracture([10, 10, 10], [1, 0, 0], (1, 1)))
        @test_throws ArgumentError embedded_fracture_domain(matrix, [f, f])
        @test_throws ArgumentError embedded_fracture_domain(matrix, f; aperture = -1)
        @test_throws ArgumentError embedded_fracture_domain(matrix, f; hydraulic_aperture = -1, permeability = 1)
        @test_throws ArgumentError embedded_fracture_domain(matrix, f; tolerance = 0)
        @test_throws ArgumentError setup_embedded_fracture_model(matrix, domain, SinglePhaseSystem(); representation = :bad)
        @test_throws ArgumentError setup_embedded_fracture_model(matrix, domain, SinglePhaseSystem(); subtract_fracture_volume = true, representation = :bad)
        zero_perm = embedded_fracture_domain(matrix, f; permeability = 0.0)
        @test all(==(0), zero_perm[:transmissibilities])
        @test all(==(0), zero_perm[:fracture_matrix_transmissibility])
        plotted = Jutul.triangulate_mesh(physical_representation(pd))
        @test size(plotted.points, 2) == 3
        @test length(plotted.mapper.Cells(ones(number_of_cells(pd)))) == size(plotted.points, 1)
    end

    @testset "Face coincidence and tensor permeability" begin
        mesh = reservoir_mesh((2, 2, 2), (2.0, 2.0, 2.0))
        K = repeat([2.0, 3.0, 4.0, 0.1, 0.2, 0.3], 1, 8)
        matrix = reservoir_domain(mesh; permeability = K)
        original_volumes = copy(matrix[:volumes])
        f = RectangleFracture([1, 1, 1], [1, 0, 0], (2, 2); permeability = 10.0, aperture = 0.01)
        domain = embedded_fracture_domain(matrix, f)
        @test number_of_cells(domain) == 4
        @test length(domain[:matrix_cells]) == 8
        @test domain[:matrix_distance] ≈ fill(0.5, 8)
        @test domain[:fracture_matrix_transmissibility] ≈ fill(1/(0.5/2 + 0.005/10), 8)
        combined = combined_fracture_domain(matrix, domain)
        @test number_of_cells(combined) == 12
        @test combined[:permeability][:, 1:8] == K
        @test all(combined[:permeability][1:3, 9:12] .== 10.0)
        @test all(combined[:permeability][4:6, 9:12] .== 0.0)
        @test combined[:neighbors][:, 1:number_of_faces(matrix)] == matrix[:neighbors]
        plotted = Jutul.triangulate_mesh(physical_representation(combined))
        @test maximum(plotted.mapper.indices.Cells) == 12
        @test length(plotted.mapper.Cells(ones(12))) == size(plotted.points, 1)
        raw = Jutul.triangulate_mesh(physical_representation(combined); flatten = false)
        @test maximum(maximum, raw.triangulation) == length(raw.points)
        adjusted = combined_fracture_domain(matrix, domain; subtract_fracture_volume = true)
        @test sum(adjusted[:volumes]) ≈ sum(matrix[:volumes])
        @test matrix[:volumes] == original_volumes
        partial = RectangleFracture([1, 1, 1], [1, 0, 0], (1.3, 0.7))
        @test sum(embedded_fracture_domain(matrix, partial)[:volumes]) ≈ 1.3*0.7*partial.properties.aperture
        tangent = PolygonFracture([[0, 0, 0], [0, 2, 0], [0, 2, 2], [0, 0, 2]])
        bd = embedded_fracture_domain(matrix, tangent)
        @test sum(bd[:volumes]) ≈ 4tangent.properties.aperture
        @test length(bd[:matrix_cells]) == 4
        @test physical_representation(bd).matrix_areas ≈ ones(4)
    end

    @testset "Unified three-dimensional mesh" begin
        mesh = reservoir_mesh((3, 3, 3), (3.0, 3.0, 3.0))
        Jutul.set_mesh_entity_tag!(mesh, Faces(), :faults, :test_fault, [1])
        matrix = reservoir_domain(mesh)
        objects = AbstractEmbeddedFracture[
            RectangleFracture([1, 1.5, 1.5], [1, 0, 0], (3, 3)),
            EllipticalFracture([1.5, 1.5, 1.5], [1, 1, 1], (0.6, 0.8); segments = 16),
            PolygonFracture([[2.5, 0, 0], [2.5, 2, 0], [2.5, 2, 1],
                [2.5, 1, 1], [2.5, 1, 2], [2.5, 0, 2]])]
        fractures = embedded_fracture_domain(matrix, objects; aperture = [0.01, 0.02, 0.03])
        combined = combined_fracture_domain(matrix, fractures)
        embedded = physical_representation(combined)
        before = deepcopy(mesh)
        converted = UnstructuredMesh(embedded)
        nm, nf, nn = number_of_cells(mesh), number_of_faces(mesh), Jutul.number_of_nodes(mesh)
        @test converted isa UnstructuredMesh{3}
        @test number_of_cells(converted) == number_of_cells(combined)
        @test number_of_faces(converted) == number_of_faces(combined)
        @test get_neighborship(converted) == combined[:neighbors]
        @test converted.node_points[1:nn] == mesh.node_points
        @test converted.faces.faces_to_nodes.vals == mesh.faces.faces_to_nodes.vals
        @test converted.boundary_faces.faces_to_nodes.vals[1:length(mesh.boundary_faces.faces_to_nodes.vals)] == mesh.boundary_faces.faces_to_nodes.vals
        @test all(f -> isempty(converted.faces.faces_to_nodes[f]), nf+1:number_of_faces(converted))
        @test Jutul.get_mesh_entity_tag(converted, Faces(), :faults, :test_fault) == [1]
        @test Jutul.grid_dims_ijk(converted) == Jutul.grid_dims_ijk(mesh)
        @test [cell_ijk(converted, c) for c in 1:number_of_cells(converted)] ==
            [cell_ijk(embedded, c) for c in 1:number_of_cells(embedded)]
        for c in 1:number_of_cells(converted)
            centroid, volume = Jutul.compute_centroid_and_measure(converted, Cells(), c)
            @test volume ≈ combined[:volumes][c]
            @test centroid ≈ combined[:cell_centroids][:, c]
            @test all(>(0), Jutul.cell_dims(converted, c))
        end
        # Each fracture boundary is closed, including partially shared edges
        # between the convex parts of a concave polygon.
        for c in nm+1:number_of_cells(converted)
            edges = Dict{Tuple{Int, Int}, Int}()
            for face in converted.boundary_faces.cells_to_faces[c]
                nodes = converted.boundary_faces.faces_to_nodes[face]
                for i in eachindex(nodes)
                    edge = (nodes[i], nodes[mod1(i+1, length(nodes))])
                    edges[edge] = get(edges, edge, 0)+1
                end
            end
            @test all(edge -> edges[edge] == get(edges, reverse(edge), 0), keys(edges))
        end
        outer = Jutul.triangulate_mesh(converted; outer = true)
        @test all(isfinite, outer.points)
        @test all(c -> c in outer.mapper.indices.Cells, nm+1:number_of_cells(converted))
        @test outer.mapper.Cells(collect(1:number_of_cells(converted))) == outer.mapper.indices.Cells
        @test !isempty(Jutul.mesh_linesegments(embedded))
        @test mesh.node_points == before.node_points
        @test mesh.faces.neighbors == before.faces.neighbors
        @test mesh.faces.cells_to_faces.vals == before.faces.cells_to_faces.vals
        @test mesh.boundary_faces.faces_to_nodes.vals == before.boundary_faces.faces_to_nodes.vals
    end

    @testset "Matrix properties and face-coincident flow paths" begin
        matrix = reservoir_domain(reservoir_mesh((2, 2, 2), (2.0, 2.0, 2.0));
            net_to_gross = 0.4, transmissibility_multiplier = fill(2.0, 12))
        matrix[:satnum, Cells()] = collect(1:8)
        matrix[:transmissibility_override, Faces()] = vcat([1e-12], fill(NaN, 11))
        before = deepcopy(matrix.data)
        object = RectangleFracture([1, 1, 1], [1, 0, 0], (2, 2))
        fractures = embedded_fracture_domain(matrix, object)
        combined = combined_fracture_domain(matrix, fractures)
        covered = findall(>(0), fractures[:matrix_face_coverage])
        @test length(covered) == 4
        @test fractures[:matrix_face_coverage][covered] ≈ ones(4)
        @test combined[:transmissibilities][covered] ≈ zeros(4) atol = 1e-25
        @test combined[:net_to_gross][9:end] == ones(4)
        parents = [first(p.matrix_cells) for p in physical_representation(fractures).patches]
        @test combined[:satnum][9:end] == matrix[:satnum][parents]
        partial = embedded_fracture_domain(matrix, RectangleFracture([1, 1, 1], [1, 0, 0], (1, 1)))
        @test sum(partial[:matrix_face_coverage]) ≈ 1
        system = SinglePhaseSystem()
        model, parameters = setup_fractured_reservoir_model(matrix, fractures, system;
            subtract_fracture_volume = true, extra_out = true)
        @test sum(reservoir_domain(model)[:volumes]) + sum(fractures[:volumes]) ≈ sum(matrix[:volumes])
        @test parameters[:Reservoir][:Transmissibilities][covered] ≈ zeros(4) atol = 1e-25
        uncovered = setdiff(1:12, covered)
        original = setup_parameters(setup_reservoir_model(matrix, system))[:Reservoir][:Transmissibilities]
        @test parameters[:Reservoir][:Transmissibilities][uncovered] == original[uncovered]
        @test Jutul.grid_dims_ijk(physical_representation(combined)) == (2, 2, 2)
        @test all(k -> isequal(matrix.data[k], before[k]), keys(before))
        base = reservoir_mesh((3, 3, 2), (3.0, 3.0, 2.0))
        nnc = setup_nnc_connections(base, [1], [18], [1e-12], [0.2])
        with_nnc = reservoir_domain(base; nnc = nnc)
        f_nnc = embedded_fracture_domain(with_nnc, RectangleFracture([1.5, 1.5, 1], [1, 0, 0], (3, 2)))
        merged = combined_fracture_domain(with_nnc, f_nnc)
        nf = number_of_faces(with_nnc)
        expected = setup_parameters(setup_reservoir_model(with_nnc, system))[:Reservoir][:Transmissibilities]
        @test merged[:transmissibilities][1:nf] == expected
        @test merged[:transmissibilities][nf] == 1e-12
        inactive = reservoir_domain(reservoir_mesh((6, 1, 1), (6.0, 1.0, 1.0); actnum = [1, 1, 0, 0, 1, 1]))
        isolated = embedded_fracture_domain(inactive, RectangleFracture([3, 0.5, 0.5], [0, 0, 1], (6, 1)))
        @test number_of_cells(isolated) == 4
        @test number_of_faces(isolated) == 2
        @test all(f -> abs(isolated[:neighbors][1, f]-isolated[:neighbors][2, f]) == 1, 1:2)
        adjacent = reservoir_domain(reservoir_mesh((2, 1, 1), (2.0, 1.0, 1.0)))
        # Fractures can meet along the shared face while their positive-area
        # patches occupy different matrix cells.
        left = PolygonFracture([[0, 0, 0.5], [1, 0, 0.5], [1, 1, 0.5], [0, 1, 0.5]])
        right = PolygonFracture([[1, 0, 0.5], [1.5, 0, 0], [1.5, 1, 0], [1, 1, 0.5]])
        meeting = embedded_fracture_domain(adjacent, [left, right])
        @test number_of_cells(meeting) == 2
        @test number_of_faces(meeting) == 1
        @test physical_representation(meeting).intersections[1].length ≈ 1
    end

    @testset "DFM network geometry and transmissibilities" begin
        for (cross, three) in ((false, false), (true, false), (true, true))
            d = efm_dfm_geometry(; cross = cross, three = three)
            @test sum(d.efm[:volumes]) ≈ sum(d.dfm[:volumes])
            @test number_of_cells(d.efm) == number_of_cells(d.dfm)
            mapping = [argmin([norm(d.efm[:cell_centroids][:, i]-d.dfm[:cell_centroids][:, j]) for j in 1:number_of_cells(d.dfm)]) for i in 1:number_of_cells(d.efm)]
            @test length(unique(mapping)) == length(mapping)
            @test d.efm[:cell_centroids] ≈ d.dfm[:cell_centroids][:, mapping]
            Tdfm = JutulDarcy.reservoir_conductivity_dfm(d.dfm, d.dfm[:permeability])
            @test efm_connection_matrix(d.efm, d.efm[:transmissibilities], mapping) ≈ efm_connection_matrix(d.dfm, Tdfm) rtol = 1e-9
            Tm = JutulDarcy.matrix_fracture_connection_conductivity(d.dfm)
            @test sort(d.efm[:fracture_matrix_transmissibility]) ≈ sort(Tm) rtol = 1e-9
        end
    end

    @testset "Coupled and single models: multiphase and geothermal" begin
        for thermal in (false, true)
            d = efm_dfm_geometry(; cross = true)
            system = ImmiscibleSystem((AqueousPhase(), LiquidPhase()); reference_densities = (1000.0, 1000.0))
            if thermal
                system = :geothermal
            end
            multi = efm_validation_case(d.matrix, d.efm, system; thermal = thermal)
            single = efm_validation_case(d.matrix, d.efm, system; thermal = thermal, representation = :single)
            volumes = copy(d.matrix[:volumes])
            solver = GenericKrylov(preconditioner = ILUZeroPreconditioner(), relative_tolerance = 1e-10)
            simarg = (info_level = -1, max_nonlinear_iterations = 20, linear_solver = solver,
                tol_cnv = 1e-10, tol_cnve = 1e-10, timesteps = :none)
            rm = simulate_reservoir(multi; simarg...)
            rs = simulate_reservoir(single; simarg...)
            @test length(rm.result.states) == length(multi.dt)
            @test length(rs.result.states) == length(single.dt)
            for (m, s) in zip(rm.result.states, rs.result.states)
                @test vcat(m[:Reservoir][:Pressure], m[:Fractures][:Pressure]) ≈ s[:Reservoir][:Pressure] rtol = 1e-7
                if thermal
                    @test vcat(m[:Reservoir][:Temperature], m[:Fractures][:Temperature]) ≈ s[:Reservoir][:Temperature] rtol = 1e-7
                else
                    @test hcat(m[:Reservoir][:Saturations], m[:Fractures][:Saturations]) ≈ s[:Reservoir][:Saturations] atol = 1e-7
                end
            end
            @test d.matrix[:volumes] == volumes
            @test physical_representation(reservoir_domain(single)) isa EFM.EmbeddedReservoirMesh
            @test !haskey(single.model.models, :Fractures)
            @test haskey(multi.model.models, :Fractures)
        end
    end

    @testset "DFM comparison on the fractured-cube validation case" begin
        coarse_temperature_error = Inf
        for (thermal, n) in ((false, (5, 5, 3)), (false, (11, 11, 7)),
                (false, (10, 10, 6)), (true, (5, 5, 3)), (true, (11, 11, 7)))
            pv_fraction = 0.5
            if thermal
                pv_fraction = 5.0
            end
            comparison = efm_dfm_comparison(n; thermal = thermal, pv_fraction = pv_fraction)
            @test length(comparison.efm.result.states) == 20
            @test length(comparison.dfm.result.states) == 20
            pressure_error, transport_error = comparison.pressure_error, comparison.transport_error
            @info "EFM/DFM fractured-cube comparison" thermal n pressure_error transport_error
            if thermal
                @test pressure_error < 0.005
                if n == (5, 5, 3)
                    @test transport_error < 2.0
                    coarse_temperature_error = transport_error
                else
                    @test transport_error < 1.0
                    @test transport_error < 0.6coarse_temperature_error
                end
            elseif all(iseven, n)
                @test pressure_error < 1e-6
                @test transport_error < 1e-6
            else
                @test pressure_error < 0.002
                # EFM retains cells straddling fractures. Transport agreement
                # on this nonconforming case is a discretization comparison,
                # and the maximum transient error need not decrease monotonically.
                @test transport_error < 0.1
            end
        end
    end
end
