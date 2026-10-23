using Jutul, JutulDarcy, Test, LinearAlgebra

function cprw_multisegment_case(; merged = false, split_wells = false)
    bar = si_unit(:bar)
    domain = reservoir_domain(CartesianMesh((3, 1, 3), (300.0, 100.0, 45.0));
        permeability = 1e-13, porosity = 0.2)
    wells = Any[setup_vertical_well(domain, i, 1; name, simple_well = false)
        for (i, name) in enumerate((:Producer, :Injector))]
    push!(wells, setup_well(domain, [3]; name = :Simple, simple_well = true))
    system = ImmiscibleSystem((LiquidPhase(), VaporPhase());
        reference_densities = [1000.0, 100.0])
    model = setup_reservoir_model(domain, system; wells, split_wells, backend = :csr)
    state0 = setup_reservoir_state(model; Pressure = 200bar, Saturations = [0.8, 0.2])
    controls = Dict(
        :Producer => ProducerControl(BottomHolePressureTarget(190bar)),
        :Injector => InjectorControl(TotalRateTarget(1e-4), [1.0, 0.0]; density = 1000.0),
        :Simple => DisabledControl(),
    )
    next_controls = copy(controls)
    next_controls[:Producer] = ProducerControl(TotalRateTarget(-1e-4))
    next_controls[:Injector] = InjectorControl(BottomHolePressureTarget(210bar), [1.0, 0.0]; density = 1000.0)
    forces = [setup_reservoir_forces(model; control) for control in (controls, next_controls)]
    case = JutulCase(model, fill(0.05*si_unit(:day), 2), forces; state0)
    return merged ? JutulDarcy.merge_similar_wells(case) : case
end

function check_multisegment_coarse_system(result)
    sim = result.extra[:simulator]
    cpr = result.extra[:config][:linear_solver].preconditioner
    map = cpr.storage.well_reservoir_map
    ncell = number_of_cells(reservoir_model(sim.model).domain)
    @test map.np == ncell + 3
    pressure = Matrix(JutulDarcy.cpr_host_matrix(cpr.storage.A_p))
    weights = Array(cpr.storage.w_p)
    ncomp = cpr.storage.number_of_components
    bz = cpr.storage.block_size
    lsys = sim.storage.LinearizedSystem
    for (key, range) in zip(map.wells, map.well_ranges)
        model = sim.model.models[key]
        well = physical_representation(model.domain)
        well isa JutulDarcy.MultiSegmentWell || continue
        group = Jutul.group_index(sim.model, key)
        offset = Jutul.group_linearized_system_offset(sim.model, key)
        # Build the cell-variable offset from the actual state field sizes;
        # segment face fluxes may precede pressure in the primary-variable list.
        variable_offset = offset
        for (name, variable) in model.primary_variables
            Jutul.associated_entity(variable) == Cells() && break
            variable_offset += length(sim.storage[key].state[name])
        end
        equation_offset = Jutul.group_linearized_system_offset(sim.model, key, Jutul.number_of_equations) +
            Jutul.get_equation_offset(model, :mass_conservation)
        layout = Jutul.matrix_layout(model.context)
        nnodes = number_of_cells(well)
        A12 = Matrix(JutulDarcy.cpr_host_matrix(lsys[1, group].jac))
        A21 = Matrix(JutulDarcy.cpr_host_matrix(lsys[group, 1].jac))
        A22 = Matrix(JutulDarcy.cpr_host_matrix(lsys[group, group].jac))
        for (i, coarse) in enumerate(range)
            if isnothing(well.multiwell)
                name = key
                nodes = 1:nnodes
            else
                name = well.multiwell.names[i]
                nodes = well.multiwell.nodes.pos[i]:(well.multiwell.nodes.pos[i + 1] - 1)
            end
            row = ncell + coarse
            cell_index(node, component) = Jutul.alignment_linear_index(node, component, nnodes, bz, layout)
            pcols = [variable_offset + cell_index(node, 1) for node in nodes]
            # Construct R_w and P_b independently of CPRW's sparse maps.
            restriction = zeros(size(A21, 1))
            for node in nodes, component in 1:ncomp
                restriction[equation_offset + cell_index(node, component)] = weights[component, row]
            end
            @test weights[:, row] ≈ vec(sum(weights[:, map.well_cells[coarse]]; dims = 2))/length(map.well_cells[coarse])
            expected = vec(sum(weights[:, 1:ncell].*reshape(vec(sum(A12[:, pcols]; dims = 2)), bz, ncell); dims = 1))
            @test pressure[1:ncell, row] ≈ expected
            if name == :Injector
                @test pressure[row, row] == 1
                @test all(iszero, pressure[row, 1:ncell])
            else
                @test pressure[row, 1:ncell] ≈ vec(transpose(restriction)*A21[:, 1:bz:end])
                @test pressure[row, row] ≈ sum(transpose(restriction)*A22[:, pcols])
            end
        end
    end
end

@testset "CPRW multisegment pressure enrichment" begin
    options = (info_level = -1, timesteps = :none, error_on_incomplete = true,
        failure_cuts_timestep = false, tol_cnv = 1e-7, tol_mb = 1e-9, tol_cnv_well = 1e-7)
    for (merged, split_wells) in ((false, false), (false, true), (true, false))
        @testset "merged=$merged, split_wells=$split_wells" begin
            case = cprw_multisegment_case(; merged, split_wells)
            reference = simulate_reservoir(case; precond = :ilu0, options...)
            result = simulate_reservoir(case; precond = :cprw, options...)
            @test length(result.states) == length(case.dt)
            for (state, ref) in zip(result.states, reference.states)
                @test state[:Pressure] ≈ ref[:Pressure] rtol = 1e-7
                @test state[:Saturations] ≈ ref[:Saturations] rtol = 1e-7
            end
            check_multisegment_coarse_system(result)
        end
    end
end
