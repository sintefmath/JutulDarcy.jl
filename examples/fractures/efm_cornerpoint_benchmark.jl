# # EFM geometry and flow on a full corner-point grid
# <tags: Validation, Discretizations>
# Usage: include this file, then call benchmark_embedded_fractures with a
# GRDECL path and optionally ACTNUM/PORO/PERM property paths. No grid data is
# downloaded or modified by this script.
using Jutul, JutulDarcy, GeoEnergyIO, Random, Test, LinearAlgebra

function efm_benchmark_matrix(grid_path;
        actnum_path = missing,
        property_paths = String[]
    )
    grid = parse_grdecl_file(grid_path; actnum_path = actnum_path, extra_paths = property_paths)
    mesh = reservoir_mesh(grid)
    cells = mesh.cell_map
    porosity = fill(0.2, number_of_cells(mesh))
    if haskey(grid, "PORO")
        porosity = vec(grid["PORO"])[cells]
    end
    permeability = fill(0.1si_unit(:darcy), 3, number_of_cells(mesh))
    if haskey(grid, "PERMX")
        permeability[1, :] = vec(grid["PERMX"])[cells]
        permeability[2, :] = permeability[1, :]
        permeability[3, :] = 0.1permeability[1, :]
        for (axis, key) in enumerate(("PERMX", "PERMY", "PERMZ"))
            if haskey(grid, key)
                permeability[axis, :] = vec(grid[key])[cells]
            end
        end
    end
    # Keep the full active geometry, including thin cells. Small positive
    # porosity only affects input cells whose PORO was zero despite ACTNUM=1.
    porosity = max.(porosity, 1e-6)
    return reservoir_domain(mesh; permeability = permeability, porosity = porosity)
end

function efm_benchmark_fractures(matrix; count = 1200, seed = 1984)
    rng = MersenneTwister(seed)
    mesh = physical_representation(matrix)
    cc = matrix[:cell_centroids]
    fractures = AbstractEmbeddedFracture[]
    for i in 1:count
        cell = rand(rng, 1:number_of_cells(matrix))
        center = cc[:, cell]
        normal = [cos(2pi*rand(rng)), sin(2pi*rand(rng)), 0.15randn(rng)]
        normal /= norm(normal)
        # Field-scale surfaces: several horizontal cell widths and multiple
        # layers, with random orientation and intersections.
        dims = Jutul.cell_dims(mesh, cell)
        width = max(100.0, 3max(dims[1], dims[2]))*(0.75+rand(rng))
        height = max(20.0, 4dims[3])*(0.75+rand(rng))
        properties = (aperture = 1e-3, permeability = 1000si_unit(:darcy), porosity = 0.5)
        if i%3 == 0
            push!(fractures, EllipticalFracture(center, normal, (height/2, width/2); direction = [0, 0, 1], segments = 32, properties...))
        elseif i%3 == 1
            push!(fractures, RectangleFracture(center, normal, (height, width); direction = [0, 0, 1], properties...))
        else
            rectangle = RectangleFracture(center, normal, (height, width); direction = [0, 0, 1], properties...)
            pts = fracture_polygon(rectangle)
            # A concave polygon with a notch tests the general polygon path.
            insert!(pts, 3, 0.6pts[2]+0.4center)
            push!(fractures, PolygonFracture(pts; properties...))
        end
    end
    return fractures
end

function benchmark_embedded_fractures(matrix::DataDomain;
        count = 1200,
        seed = 1984,
        representations = (:single, :coupled),
        solve = true,
        amg = :hypre,
        info_level = -1,
        scaling = :none
    )
    fractures = efm_benchmark_fractures(matrix; count = count, seed = seed)
    mesh = physical_representation(matrix)
    nodes = copy(mesh.node_points)
    neighbors = copy(matrix[:neighbors])
    volumes = copy(matrix[:volumes])
    # Warm up geometry methods independently of the measured large workload.
    embedded_fracture_domain(matrix, fractures[1:3])
    GC.gc()
    timed = @timed embedded_fracture_domain(matrix, fractures)
    domain = timed.value
    fmesh = physical_representation(domain)
    @test mesh.node_points == nodes
    @test matrix[:neighbors] == neighbors
    @test matrix[:volumes] == volumes
    @test all(isfinite, domain[:transmissibilities])
    @test all(>=(0), domain[:transmissibilities])
    @test all(isfinite, domain[:fracture_matrix_transmissibility])
    @test all(>(0), domain[:volumes])
    @test length(unique(domain[:fracture_id])) == count
    report = (matrix_cells = number_of_cells(matrix), fractures = count,
        fracture_cells = number_of_cells(domain), fracture_connections = number_of_faces(domain),
        matrix_connections = length(fmesh.geometry.matrix_cells), intersection_segments = length(fmesh.geometry.intersections),
        preprocessing_seconds = timed.time, preprocessing_allocated_bytes = timed.bytes)
    println(report)
    flow = Dict{Symbol, Any}()
    if solve
        cc = matrix[:cell_centroids]
        lo, hi = argmin(cc[1, :]), argmax(cc[1, :])
        bc = [FlowBoundaryCondition(matrix, lo, 110si_unit(:bar)),
              FlowBoundaryCondition(matrix, hi, 90si_unit(:bar))]
        system = SinglePhaseSystem(AqueousPhase(); reference_density = 1000.0)
        for representation in representations
            setup = @timed setup_embedded_fracture_model(matrix, domain, system;
                representation = representation, block_backend = false, backend = :csc)
            model = setup.value
            state0 = setup_reservoir_state(model; Pressure = 100si_unit(:bar))
            forces = setup_reservoir_forces(model; bc = bc)
            case = JutulCase(model, fill(si_unit(:day), 3), forces; state0 = state0)
            if amg == :hypre
                preconditioner = BoomerAMGPreconditioner(; AggNumLevels = 0, RelaxType = 6)
            elseif amg == :jutul
                preconditioner = AMGPreconditioner(:hmis; coarse_size = 32)
            else
                throw(ArgumentError("amg must be :hypre or :jutul"))
            end
            solver = GenericKrylov(:bicgstab; preconditioner = preconditioner, scaling = scaling, relative_tolerance = 1e-8, max_iterations = 300)
            println((representation = representation, status = :solving, cells = number_of_cells(reservoir_domain(model))))
            flush(stdout)
            result = @timed simulate_reservoir(case; linear_solver = solver, info_level = info_level,
                max_nonlinear_iterations = 12, max_timestep_cuts = 5, failure_cuts_timestep = false)
            states = result.value.result.states
            @test length(states) == 3
            @test all(isfinite, states[end][:Reservoir][:Pressure])
            @test maximum(abs.(states[end][:Reservoir][:Pressure] .- state0[:Reservoir][:Pressure])) > 1.0
            stats = Jutul.report_stats(result.value.result.reports)
            flow[representation] = (setup_seconds = setup.time, solve_seconds = result.time,
                linear_iterations = stats.linear_iterations, nonlinear_iterations = stats.newtons,
                result = result.value)
            println((representation = representation, setup_seconds = setup.time,
                solve_seconds = result.time, linear_iterations = stats.linear_iterations,
                nonlinear_iterations = stats.newtons))
        end
        if haskey(flow, :single) && haskey(flow, :coupled)
            single = flow[:single].result.result.states[end][:Reservoir][:Pressure]
            coupled = flow[:coupled].result.result.states[end]
            p = vcat(coupled[:Reservoir][:Pressure], coupled[:Fractures][:Pressure])
            @test maximum(abs.(single-p))/maximum(abs, p) < 1e-5
        end
    end
    return (; report, flow, domain)
end

function benchmark_embedded_fractures(grid_path::AbstractString; actnum_path = missing, property_paths = String[], kwarg...)
    matrix = efm_benchmark_matrix(grid_path; actnum_path = actnum_path, property_paths = property_paths)
    return benchmark_embedded_fractures(matrix; kwarg...)
end
