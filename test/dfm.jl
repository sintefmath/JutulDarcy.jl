using Jutul, JutulDarcy, Test, LinearAlgebra
import Jutul.CutCellMeshes: PlaneCut, cut_mesh

Darcy = si_unit(:darcy)
bar   = si_unit(:bar)
meter = si_unit(:meter)
kg    = si_unit(:kilogram)
year  = si_unit(:year)

## ── Helper: build a simple fractured domain ─────────────────────────────
#  Returns (matrix, fractures, fracture_mesh, matrix_mesh) for a unit cube
#  with a single fracture at x = L/2.
function setup_simple_dfm(;
        n = (3, 3, 3),
        L = (1.0, 1.0, 1.0),
        matrix_permeability = 1e-2Darcy,
        matrix_porosity     = 0.1,
        fracture_aperture   = 1e-3meter,
        fracture_porosity   = 0.5,
        intersection_strategy = :star_delta,
    )
    dims = n isa Integer ? (n, n, n) : Tuple(n)
    Lxyz = L isa Number ? (L, L, L) : Tuple(L)

    mesh = UnstructuredMesh(CartesianMesh(dims, Lxyz))
    ijk0 = reinterpret(reshape, Int,
        map(c -> cell_ijk(mesh, c), 1:number_of_cells(mesh)))

    cuts = [Jutul.PlaneCut(
        [Lxyz[1] / 2, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    )]
    cut, info = Jutul.cut_mesh(mesh, cuts; extra_out = true)
    ijk = ijk0[:, info[:cell_index]]

    fracture_faces = findall(info[:face_index] .== 0)
    @assert !isempty(fracture_faces) "No fracture faces found after cut"

    fracture_mesh = Jutul.EmbeddedMesh(cut, fracture_faces;
        intersection_strategy = intersection_strategy)

    matrix = reservoir_domain(cut;
        permeability = matrix_permeability,
        porosity     = matrix_porosity,
    )

    fractures = JutulDarcy.fracture_domain(fracture_mesh, matrix;
        aperture = fracture_aperture,
        porosity = fracture_porosity,
    )

    return (; matrix, fractures, fracture_mesh, cut, ijk)
end

## ── Helper: build a cross-shaped fracture network (fractures at x=L/2 and y=L/2) ─
function setup_cross_dfm(;
        n = (3, 3, 3),
        L = (1.0, 1.0, 1.0),
        matrix_permeability = 1e-2Darcy,
        matrix_porosity     = 0.1,
        fracture_aperture   = 1e-3meter,
        fracture_porosity   = 0.5,
        intersection_strategy = :star_delta,
    )
    dims = n isa Integer ? (n, n, n) : Tuple(n)
    Lxyz = L isa Number ? (L, L, L) : Tuple(L)

    mesh = UnstructuredMesh(CartesianMesh(dims, Lxyz))

    cuts = [
        Jutul.PlaneCut([Lxyz[1] / 2, 0.0, 0.0], [1.0, 0.0, 0.0]),
        Jutul.PlaneCut([0.0, Lxyz[2] / 2, 0.0], [0.0, 1.0, 0.0]),
    ]
    cut, info = Jutul.cut_mesh(mesh, cuts; extra_out = true)
    fracture_faces = findall(info[:face_index] .== 0)

    fracture_mesh = Jutul.EmbeddedMesh(cut, fracture_faces;
        intersection_strategy = intersection_strategy)

    matrix = reservoir_domain(cut;
        permeability = matrix_permeability,
        porosity     = matrix_porosity,
    )
    fractures = JutulDarcy.fracture_domain(fracture_mesh, matrix;
        aperture = fracture_aperture,
        porosity = fracture_porosity,
    )

    return (; matrix, fractures, fracture_mesh, cut)
end

## ── Helper: build FD + DFM case pair for validation ─────────────────────
#  Creates a fractured cube with a single fracture at x = L/2 and wells in
#  opposing corners. Returns ready-to-simulate JutulCase objects for both
#  full-dimensional (FD) and discrete fracture model (DFM) approaches.
function setup_fractured_cube_test(;
        n = (5, 5, 3),
        L = (10.0, 10.0, 3.0),
        matrix_permeability = 1e-2Darcy,
        matrix_porosity = 0.05,
        fracture_aperture = 1e-4,
        fracture_porosity = 0.5,
        pv_frac = 0.5,
        nstep = 20,
        total_time = 0.1year,
        fluid_model = :two_phase,
        temperature = convert_to_si(90.0, :Celsius),
        injection_temperature = convert_to_si(10.0, :Celsius),
    )
    dims0 = n isa Integer ? (n, n, n) : Tuple(n)
    Lxyz  = L isa Number  ? (L, L, L) : Tuple(L)

    ## ── Full-dimensional (FD) mesh ──────────────────────────────────
    fd_axis_widths = ntuple(3) do d
        nd = dims0[d]; Ld = Lxyz[d]
        xc = collect(range(0.0, Ld; length = nd + 1))
        xf = Ld / 2
        x = sort(vcat(xc, xf .- fracture_aperture / 2, xf .+ fracture_aperture / 2))
        x = unique(filter(xi -> 0.0 <= xi <= Ld, x))
        diff(x)
    end
    fd_cart_dims = map(length, fd_axis_widths)

    fd_mesh = UnstructuredMesh(CartesianMesh(fd_cart_dims, fd_axis_widths))
    fd_ijk = reinterpret(reshape, Int,
        map(c -> cell_ijk(fd_mesh, c), 1:number_of_cells(fd_mesh)))

    fd_matrix = reservoir_domain(fd_mesh;
        permeability = matrix_permeability, porosity = matrix_porosity)
    is_frac = [any(isapprox.(cell_dims(fd_mesh, c), fracture_aperture; rtol = 1e-2))
               for c in 1:number_of_cells(fd_mesh)]
    fd_matrix[:permeability][is_frac] .= fracture_aperture^2 / 12
    fd_matrix[:porosity][is_frac] .= fracture_porosity

    radius = 0.1
    fd_inj = setup_well(fd_matrix,
        findall(fd_ijk[1, :] .== 1 .&& fd_ijk[2, :] .== 1);
        name = :Injector, simple_well = false, dir = :z, radius = radius)
    fd_prod = setup_well(fd_matrix,
        findall(fd_ijk[1, :] .== maximum(fd_ijk[1, :]) .&& fd_ijk[2, :] .== maximum(fd_ijk[2, :]));
        name = :Producer, simple_well = false, dir = :z, radius = radius)
    fd_pv = sum(pore_volume(fd_matrix))

    ## ── Lower-dimensional (DFM) mesh ────────────────────────────────
    ld_mesh = UnstructuredMesh(CartesianMesh(dims0, Lxyz))
    ld_ijk = reinterpret(reshape, Int,
        map(c -> cell_ijk(ld_mesh, c), 1:number_of_cells(ld_mesh)))

    cuts = [
        PlaneCut([Lxyz[1] / 2, 0.0, 0.0], [1.0, 0.0, 0.0])
        PlaneCut([0.0, Lxyz[2] / 2, 0.0], [0.0, 1.0, 0.0])
        PlaneCut([0.0, 0.0, Lxyz[3] / 2], [0.0, 0.0, 1.0])
        ]
    ld_cut_mesh, info = cut_mesh(ld_mesh, cuts; extra_out = true)
    fracture_faces = findall(info[:face_index] .== 0)
    fracture_mesh = Jutul.EmbeddedMesh(ld_cut_mesh, fracture_faces;
        intersection_strategy = :keep)
    ld_matrix = reservoir_domain(ld_cut_mesh;
        permeability = matrix_permeability, porosity = matrix_porosity)
    ld_fractures = JutulDarcy.fracture_domain(fracture_mesh, ld_matrix;
        aperture = fracture_aperture, porosity = fracture_porosity)
    ijk = ld_ijk[:, info[:cell_index]]

    function sorted_well_cells(domain, ijk, d1_val, d2_val)
        cells = findall(ijk[1, :] .== d1_val .&& ijk[2, :] .== d2_val)
        x = domain[:cell_centroids][:, cells]
        cells[sortperm(vec(sum(x .^ 2; dims = 1)))]
    end

    inj_cells = sorted_well_cells(ld_matrix, ijk, 1, 1)
    ld_inj = setup_well(ld_matrix, inj_cells;
        name = :Injector, dir = :z, radius = radius, simple_well = false)
    prod_cells = sorted_well_cells(ld_matrix, ijk, maximum(ijk[1, :]), maximum(ijk[2, :]))
    ld_prod = setup_well(ld_matrix, prod_cells;
        name = :Producer, dir = :z, radius = radius, simple_well = false)
    ld_pv = sum(pore_volume(ld_matrix)) + sum(pore_volume(ld_fractures))

    ## ── Fluid model & cases ─────────────────────────────────────────
    thermal = false
    if fluid_model == :two_phase
        rho = (800.0kg / meter^3, 1000.0kg / meter^3)
        system = ImmiscibleSystem((AqueousPhase(), LiquidPhase());
            reference_densities = rho)
        init = Dict(:Pressure => 50bar, :Saturations => [0.0, 1.0])
        mix = [1.0, 0.0]
    elseif fluid_model == :geothermal
        rho = (1000.0kg / meter^3,)
        system = :geothermal
        init = Dict(:Pressure => 50bar, :Temperature => temperature)
        mix = [1.0]
        thermal = true
    end

    function make_case(matrix, wells, fractures, pv)
        if isnothing(fractures)
            model = setup_reservoir_model(matrix, system;
                thermal = thermal, wells = wells)
        else
            model = JutulDarcy.setup_fractured_reservoir_model(
                matrix, fractures, system;
                wells = wells, thermal = thermal)
        end
        dt = fill(total_time / nstep, nstep)
        rate = pv_frac * pv / sum(dt)
        inj_kwargs = Dict{Symbol, Any}(:density => rho[1])
        if thermal
            inj_kwargs[:temperature] = injection_temperature
        end
        ctrl_inj = InjectorControl(TotalRateTarget(rate), mix; inj_kwargs...)
        ctrl_prod = ProducerControl(BottomHolePressureTarget(10.0bar))
        control = Dict(:Injector => ctrl_inj, :Producer => ctrl_prod)
        forces = setup_reservoir_forces(model; control = control)
        state0 = setup_reservoir_state(model; init...)
        return JutulCase(model, dt, forces; state0 = state0)
    end

    case_fd  = make_case(fd_matrix, [fd_inj, fd_prod], nothing, fd_pv)
    case_dfm = make_case(ld_matrix, [ld_inj, ld_prod], ld_fractures, ld_pv)

    return (case_fd = case_fd, case_dfm = case_dfm,
            dfm_domain = ld_matrix, dfm_fractures = ld_fractures,
            fd_domain = fd_matrix)
end

## ── Mapping helpers (from dfm_validation.jl) ────────────────────────────
function map_dfm_to_fd(fd_domain, dfm_domain, dfm_fractures)
    x_fd = fd_domain[:cell_centroids]
    nearest(x) = last(findmin(vec(sum((x .- x_fd) .^ 2; dims = 1))))
    ld2fd_m = [nearest(c) for c in eachcol(dfm_domain[:cell_centroids])]
    ld2fd_f = [nearest(c) for c in eachcol(dfm_fractures[:cell_centroids])]
    return ld2fd_m, ld2fd_f
end

function reconstruct_full_state(state_m, state_f, state_fd, ld2fd_m, ld2fd_f)
    full = deepcopy(state_fd)
    nc = length(state_m[:Pressure])
    for k in keys(state_m)
        haskey(full, k) || continue
        v, vm, vf = full[k], state_m[k], state_f[k]
        if vm isa Vector && length(vm) == nc
            v[ld2fd_m] .= vm
            v[ld2fd_f] .= vf
        elseif vm isa AbstractArray && size(vm, 2) == nc
            v[:, ld2fd_m] .= vm
            v[:, ld2fd_f] .= vf
        end
        full[k] = v
    end
    return full
end

## ═══════════════════════════════════════════════════════════════════════
@testset "DFM" begin

    @testset "fracture_domain (standalone)" begin
        d = setup_simple_dfm()
        fmesh = d.fracture_mesh
        nc = number_of_cells(fmesh)

        # Basic domain - no matrix connection
        fd = JutulDarcy.fracture_domain(fmesh; aperture = 1e-3meter, porosity = 0.5)
        @test fd isa DataDomain
        @test haskey(fd, :permeability)
        @test haskey(fd, :aperture)
        @test haskey(fd, :volumes)
        @test haskey(fd, :areas)

        # Cubic law: k = a^2/12
        a = 1e-3
        expected_perm = a^2 / 12
        @test all(fd[:permeability] .≈ expected_perm)
        @test all(fd[:aperture] .≈ a)

        # Volumes should be surface area * aperture
        geo = tpfv_geometry(fmesh)
        @test all(fd[:volumes][1:nc] .> 0)
    end

    @testset "fracture_domain input validation" begin
        d = setup_simple_dfm()
        fmesh = d.fracture_mesh

        @test_throws ArgumentError JutulDarcy.fracture_domain(fmesh; aperture = NaN)
        @test_throws ArgumentError JutulDarcy.fracture_domain(fmesh; aperture = -1.0)
        @test_throws ArgumentError JutulDarcy.fracture_domain(fmesh; hydralic_aperture = Inf)
        @test_throws ArgumentError JutulDarcy.fracture_domain(fmesh; hydralic_aperture = -0.1)
    end

    @testset "fracture_domain (with matrix connection)" begin
        d = setup_simple_dfm()
        fractures = d.fractures
        fmc = JutulDarcy.FractureMatrixConnection()

        # Should have FractureMatrixConnection entity
        @test haskey(fractures.entities, fmc)
        n_conn = count_entities(fractures, fmc)
        @test n_conn > 0

        # Should have matrix connection data
        @test haskey(fractures, :connection_cells)
        @test haskey(fractures, :matrix_faces)
        @test haskey(fractures, :matrix_cells)
        @test haskey(fractures, :matrix_cell_centroids)
        @test haskey(fractures, :matrix_permeability)
        @test haskey(fractures, :matrix_porosity)
        @test haskey(fractures, :cell_normals)

        conn_cells = fractures[:connection_cells, fmc]
        matrix_faces = fractures[:matrix_faces, fmc]
        matrix_cells = fractures[:matrix_cells, fmc]

        # Connection arrays should have consistent lengths
        @test length(conn_cells) == n_conn
        @test length(matrix_faces) == n_conn
        @test length(matrix_cells) == n_conn

        # Cell normals should be unit vectors for non-intersection cells
        cn = fractures[:cell_normals, Cells()]
        for c in conn_cells
            n = cn[:, c]
            @test norm(n) ≈ 1.0
        end
    end

    @testset "fracture_domain (cross-shaped, intersections)" begin
        d = setup_cross_dfm(intersection_strategy = :star_delta)
        fractures = d.fractures
        fmesh = d.fracture_mesh

        # Should have intersection neighbors
        @test !isempty(fmesh.intersection_neighbors)

        # Connection cells should exclude intersection cells
        fmc = JutulDarcy.FractureMatrixConnection()
        conn_cells = fractures[:connection_cells, fmc]
        if hasproperty(fmesh, :intersection_cells)
            for c in conn_cells
                @test !(c in fmesh.intersection_cells)
            end
        end
    end

    @testset "matrix_fracture_connection_conductivity" begin
        d = setup_simple_dfm()
        fractures = d.fractures

        T = JutulDarcy.matrix_fracture_connection_conductivity(fractures, :permeability)
        @test length(T) > 0
        @test all(T .> 0)
        @test all(isfinite, T)

        T_rock = JutulDarcy.matrix_fracture_connection_conductivity(fractures, :rock_thermal_conductivity)
        @test length(T_rock) == length(T)
        @test all(T_rock .>= 0)
        @test all(isfinite, T_rock)

        T_fluid = JutulDarcy.matrix_fracture_connection_conductivity(fractures, :fluid_thermal_conductivity)
        @test length(T_fluid) == length(T)
        @test all(T_fluid .>= 0)
        @test all(isfinite, T_fluid)
    end

    @testset "matrix_fracture_connection_conductivity (direct)" begin
        # Simple geometry: fracture at x = 0.5, matrix cell center at x = 0.75
        aperture = 1e-3
        area = 1.0
        normal = [1.0, 0.0, 0.0]
        xf = [0.5, 0.5, 0.5]
        xm = [0.75, 0.5, 0.5]
        Kf = 1e-6
        Km = 1e-12

        T = JutulDarcy.matrix_fracture_connection_conductivity(
            aperture, area, normal, xf, xm, Kf, Km)
        @test T > 0
        @test isfinite(T)

        # Symmetry: swapping fracture/matrix permeability should change T
        T_swap = JutulDarcy.matrix_fracture_connection_conductivity(
            aperture, area, normal, xf, xm, Km, Kf)
        @test T_swap > 0
        @test T_swap != T

        # Doubling aperture should change transmissibility
        T_2a = JutulDarcy.matrix_fracture_connection_conductivity(
            2 * aperture, area, normal, xf, xm, Kf, Km)
        @test T_2a != T
    end

    @testset "setup_fractured_reservoir_model (immiscible, no wells)" begin
        d = setup_simple_dfm()
        sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))

        model = JutulDarcy.setup_fractured_reservoir_model(d.matrix, d.fractures, sys)

        @test model isa Jutul.MultiModel
        @test haskey(model.models, :Reservoir)
        @test haskey(model.models, :Fractures)

        # Should have cross-terms
        @test !isempty(model.cross_terms)

        # At least one MatrixFromFractureFlowCT
        has_flow_ct = any(ct -> ct.cross_term isa JutulDarcy.MatrixFromFractureFlowCT, model.cross_terms)
        @test has_flow_ct

        # Fracture model should use DFM transmissibilities
        frac_model = model.models[:Fractures]
        params = Jutul.get_parameters(frac_model)
        @test haskey(params, :Transmissibilities)
        @test params[:Transmissibilities] isa JutulDarcy.TransmissibilitiesDFM
    end

    @testset "setup_fractured_reservoir_model (geothermal, no wells)" begin
        d = setup_simple_dfm()
        sys = :geothermal

        model = JutulDarcy.setup_fractured_reservoir_model(d.matrix, d.fractures, sys;
            thermal = true)

        @test model isa Jutul.MultiModel
        @test haskey(model.models, :Reservoir)
        @test haskey(model.models, :Fractures)

        # Should have thermal cross-terms
        has_thermal_ct = any(ct -> ct.cross_term isa JutulDarcy.MatrixFromFractureThermalCT, model.cross_terms)
        @test has_thermal_ct

        has_flow_ct = any(ct -> ct.cross_term isa JutulDarcy.MatrixFromFractureFlowCT, model.cross_terms)
        @test has_flow_ct

        # Fracture model should use DFM transmissibilities and thermal conductivities
        frac_model = model.models[:Fractures]
        params = Jutul.get_parameters(frac_model)
        @test haskey(params, :Transmissibilities)
        @test params[:Transmissibilities] isa JutulDarcy.TransmissibilitiesDFM
        @test haskey(params, :RockThermalConductivities)
        @test params[:RockThermalConductivities] isa JutulDarcy.RockThermalConductivitiesDFM
        @test haskey(params, :FluidThermalConductivities)
        @test params[:FluidThermalConductivities] isa JutulDarcy.FluidThermalConductivitiesDFM
    end

    @testset "block_fracture_face_connections!" begin
        d = setup_simple_dfm()
        sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))

        model = JutulDarcy.setup_fractured_reservoir_model(d.matrix, d.fractures, sys)
        matrix_model = model.models[:Reservoir]
        matrix_domain = matrix_model.data_domain

        fmc = JutulDarcy.FractureMatrixConnection()
        matrix_faces = unique(d.fractures[:matrix_faces, fmc])

        # Transmissibilities on fracture faces should be zero
        T = matrix_domain[:transmissibilities, Faces()]
        for f in matrix_faces
            @test T[f] == 0.0
        end

        # Non-fracture faces should retain positive transmissibilities
        all_faces = 1:number_of_faces(matrix_domain)
        non_frac = setdiff(all_faces, matrix_faces)
        @test any(T[non_frac] .> 0)
    end

    @testset "adjust_matrix_cell_volumes!" begin
        d = setup_simple_dfm()
        sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))

        # Get volumes before and after adjustment
        original_volumes = copy(d.matrix[:volumes, Cells()])
        model = JutulDarcy.setup_fractured_reservoir_model(d.matrix, d.fractures, sys)
        adjusted_volumes = model.models[:Reservoir].data_domain[:volumes, Cells()]

        fmc = JutulDarcy.FractureMatrixConnection()
        matrix_cells = d.fractures[:matrix_cells, fmc]
        affected = unique(matrix_cells)

        # Affected cells should have reduced volumes
        for c in affected
            @test adjusted_volumes[c] < original_volumes[c]
        end

        # Unaffected cells should keep original volumes
        unaffected = setdiff(1:length(original_volumes), affected)
        for c in unaffected
            @test adjusted_volumes[c] ≈ original_volumes[c]
        end
    end

    @testset "setup_fractured_reservoir_model (with wells)" begin
        d = setup_simple_dfm(n = (3, 3, 3), L = (3.0, 3.0, 3.0))
        matrix = d.matrix
        cut = d.cut

        # ijk = reinterpret(reshape, Int,
        #     map(c -> cell_ijk(cut, c), 1:number_of_cells(cut)))
        ijk = d.ijk

        inj_cells = findall(ijk[1, :] .== 1 .&& ijk[2, :] .== 1)
        prod_cells = findall(ijk[1, :] .== maximum(ijk[1, :]) .&& ijk[2, :] .== maximum(ijk[2, :]))

        inj = setup_well(matrix, inj_cells; name = :Injector, simple_well = false)
        prod = setup_well(matrix, prod_cells; name = :Producer, simple_well = false)

        sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))
        model = JutulDarcy.setup_fractured_reservoir_model(
            d.matrix, d.fractures, sys; wells = [inj, prod])

        @test haskey(model.models, :Injector)
        @test haskey(model.models, :Producer)
        @test haskey(model.models, :Reservoir)
        @test haskey(model.models, :Fractures)
    end

    @testset "immiscible simulation runs" begin
        setup = setup_fractured_cube_test(fluid_model = :two_phase)

        sim_args = (info_level = -1, max_nonlinear_iterations = 15, max_timestep_cuts = 25)
        res_fd  = simulate_reservoir(setup.case_fd;  sim_args...)
        res_dfm = simulate_reservoir(setup.case_dfm; sim_args...)

        @test length(res_fd.result.states) > 0
        @test length(res_dfm.result.states) > 0

        # Map DFM cells to FD cells and reconstruct full state for comparison
        ld2fd_m, ld2fd_f = map_dfm_to_fd(setup.fd_domain, setup.dfm_domain, setup.dfm_fractures)

        states_fd  = map(s -> s[:Reservoir], res_fd.result.states)
        states_dfm = map(s -> s[:Reservoir], res_dfm.result.states)
        states_frac = map(s -> s[:Fractures], res_dfm.result.states)

        nstates = min(length(states_fd), length(states_dfm))
        states_recon = [
            reconstruct_full_state(states_dfm[i], states_frac[i], states_fd[i], ld2fd_m, ld2fd_f)
            for i in 1:nstates
        ]

        # Compare final-step saturation: RMS difference should be small
        S_fd  = states_fd[nstates][:Saturations]
        S_dfm = states_recon[nstates][:Saturations]
        rms_sat = sqrt(sum((S_fd .- S_dfm) .^ 2) / length(S_fd))
        @test rms_sat < 1e-3

        # Compare final-step pressure: relative RMS should be small
        p_fd  = states_fd[nstates][:Pressure]
        p_dfm = states_recon[nstates][:Pressure]
        rms_p = sqrt(sum(((p_fd .- p_dfm) ./ p_fd) .^ 2) / length(p_fd))
        @test rms_p < 1e-3
    end

    @testset "geothermal simulation runs" begin
        setup = setup_fractured_cube_test(
            fluid_model = :geothermal,
            pv_frac = 5.0,
        )

        sim_args = (info_level = -1, initial_dt = 5.0)
        res_fd  = simulate_reservoir(setup.case_fd;  sim_args...)
        res_dfm = simulate_reservoir(setup.case_dfm; sim_args...)

        @test length(res_fd.result.states) > 0
        @test length(res_dfm.result.states) > 0

        # Map DFM cells to FD cells and reconstruct full state for comparison
        ld2fd_m, ld2fd_f = map_dfm_to_fd(setup.fd_domain, setup.dfm_domain, setup.dfm_fractures)

        states_fd  = map(s -> s[:Reservoir], res_fd.result.states)
        states_dfm = map(s -> s[:Reservoir], res_dfm.result.states)
        states_frac = map(s -> s[:Fractures], res_dfm.result.states)

        nstates = min(length(states_fd), length(states_dfm))
        states_recon = [
            reconstruct_full_state(states_dfm[i], states_frac[i], states_fd[i], ld2fd_m, ld2fd_f)
            for i in 1:nstates
        ]

        # Compare final-step temperature: RMS difference should be small (in Kelvin)
        T_fd  = states_fd[nstates][:Temperature]
        T_dfm = states_recon[nstates][:Temperature]
        rms_T = sqrt(sum((T_fd .- T_dfm) .^ 2) / length(T_fd))
        @test rms_T < 7e-2 # within 10 K

        # Compare final-step pressure: relative RMS should be small
        p_fd  = states_fd[nstates][:Pressure]
        p_dfm = states_recon[nstates][:Pressure]
        rms_p = sqrt(sum(((p_fd .- p_dfm) ./ p_fd) .^ 2) / length(p_fd))
        @test rms_p < 1e-3
    end

    @testset "cross-shaped fracture network model" begin
        d = setup_cross_dfm(n = (3, 3, 3), L = (1.0, 1.0, 1.0))
        sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))

        model = JutulDarcy.setup_fractured_reservoir_model(d.matrix, d.fractures, sys)

        @test model isa Jutul.MultiModel
        @test haskey(model.models, :Fractures)

        state0 = setup_reservoir_state(model,
            Pressure = 100bar, Saturations = [0.0, 1.0])

        nstep = 5
        dt = fill(0.01year, nstep)
        forces = setup_reservoir_forces(model)

        sim = Simulator(model, state0 = state0)
        cfg = simulator_config(sim, info_level = -1)
        result = simulate(sim, dt, forces = forces, config = cfg)

        @test length(result.states) == nstep
    end

    @testset "DFM cross-term types" begin
        mc = [1, 2, 3]
        fc = [4, 5, 6]

        ct_flow = JutulDarcy.MatrixFromFractureFlowCT(mc, fc)
        @test ct_flow.matrix_cells == mc
        @test ct_flow.fracture_cells == fc
        @test Jutul.symmetry(ct_flow) isa Jutul.CTSkewSymmetry

        ct_thermal = JutulDarcy.MatrixFromFractureThermalCT(mc, fc)
        @test ct_thermal.matrix_cells == mc
        @test ct_thermal.fracture_cells == fc

        wc = [1, 2]
        fwc = [3, 4]
        ct_well = JutulDarcy.FracturesFromWellFlowCT(fwc, wc)
        @test ct_well.fracture_cells == fwc
        @test ct_well.well_cells == wc
    end

    @testset "DFM parameter types" begin
        # FractureMatrixTransmissibility
        fmt = JutulDarcy.FractureMatrixTransmissibility()
        @test Jutul.minimum_value(fmt) == 0.0
        @test Jutul.associated_entity(fmt) == JutulDarcy.FractureMatrixConnection()

        # FractureMatrixGravityDifference
        fmg = JutulDarcy.FractureMatrixGravityDifference()
        @test Jutul.minimum_value(fmg) == -Inf
        @test Jutul.associated_entity(fmg) == JutulDarcy.FractureMatrixConnection()

        # TransmissibilitiesDFM
        tdfm = JutulDarcy.TransmissibilitiesDFM()
        @test Jutul.minimum_value(tdfm) == 0.0
        @test Jutul.associated_entity(tdfm) == Faces()

        # FractureWellIndices
        fwi = JutulDarcy.FractureWellIndices()
        @test Jutul.minimum_value(fwi) == 0.0
        @test Jutul.associated_entity(fwi) == JutulDarcy.FracturePerforations()
    end

    @testset "reservoir_conductivity_dfm" begin
        d = setup_simple_dfm()
        fractures = d.fractures
        K = fractures[:permeability]

        T = JutulDarcy.reservoir_conductivity_dfm(fractures, K)
        nf = number_of_faces(fractures)
        @test length(T) == nf
        @test all(T .>= 0)
        @test all(isfinite, T)
    end

end
