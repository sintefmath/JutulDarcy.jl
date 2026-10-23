using Jutul, JutulDarcy, Test, JLArrays, CUDA, SparseArrays
using Jutul.KernelExecution: KernelAbstractions

@testset "CPRW multisegment reduction on device backends" begin
    options = (info_level = -1, timesteps = :none, error_on_incomplete = true,
        tol_cnv = 1e-7, tol_mb = 1e-9, tol_cnv_well = 1e-7)
    backends = Any[KernelAbstractions.CPU(), JLBackend()]
    if CUDA.functional()
        CUDA.allowscalar(false)
        push!(backends, CUDA.CUDABackend())
    end
    for merged in (false, true)
        case = cprw_multisegment_case(; merged)
        result = simulate_reservoir(case; precond = :cprw, options...)
        sim = result.extra[:simulator]
        cpr = result.extra[:config][:linear_solver].preconditioner
        A = JutulDarcy.reservoir_jacobian(sim.storage.LinearizedSystem)
        original_map = cpr.storage.well_reservoir_map
        # The final convergence check can assemble a newer Jacobian without
        # taking another linear solve. Compare both reductions at that Jacobian.
        JutulDarcy.update_cpr_pressure_system!(cpr.storage.A_p, cpr.pressure_precond,
            A, cpr.storage.w_p, reservoir_model(sim.model).context,
            Jutul.default_executor(), original_map, sim.model, sim.storage)
        reference = Matrix(cpr.storage.A_p)
        ncell = size(A, 1)
        fixed_rows = [row for row in (ncell + 1):original_map.np if reference[row, row] == 1]
        @test length(fixed_rows) == 1
        for backend in backends
            @testset "merged=$merged, $backend" begin
                context = KernelAbstractionsContext(backend;
                    matrix_layout = Jutul.matrix_layout(reservoir_model(sim.model).context))
                adapt(array) = JutulDarcy.Adapt.adapt(context, array)
                matrix = Jutul.StaticSparsityMatrixCSR(
                    adapt(copy(nonzeros(A))), adapt(Jutul.colvals(A)), adapt(A.rowptr),
                    size(A)..., backend; minbatch = 32)
                map = merge(original_map, (
                    nzmap_11 = Int[], nzmap_12 = Int[], nzmap_21 = Int[], nzmap_22 = Int[],
                    nzval_12 = adapt(copy(original_map.nzval_12)),
                    nzval_21 = adapt(copy(original_map.nzval_21)),
                    nzval_22 = adapt(copy(original_map.nzval_22)),
                    source_values = nothing, backend_data = Ref{Any}(nothing),
                ))
                pressure = JutulDarcy.pressure_matrix_from_global_jacobian(matrix,
                    Float64, sim.storage.LinearizedSystem, map)
                weights = adapt(copy(cpr.storage.w_p))
                for (key, range) in zip(map.wells, map.well_ranges)
                    if physical_representation(sim.model.models[key].domain) isa MultiSegmentWell
                        fill!(view(weights, :, (ncell + first(range)):(ncell + last(range))), 0)
                        JutulDarcy.cprw_average_weights!(weights, map, range, ncell, context)
                    end
                end
                @test Array(weights) ≈ cpr.storage.w_p
                for factor in (1.0, 1.25)
                    nonzeros(matrix) .*= factor
                    for values in (map.nzval_12, map.nzval_21, map.nzval_22)
                        values .*= factor
                    end
                    JutulDarcy.update_cpr_pressure_system!(pressure, cpr.pressure_precond,
                        matrix, weights, context, Jutul.default_executor(), map, sim.model, sim.storage)
                    expected = factor.*reference
                    for row in fixed_rows
                        expected[row, row] = 1
                    end
                    @test Matrix(JutulDarcy.cpr_host_matrix(pressure)) ≈ expected rtol = 1e-12
                end
            end
        end
        policies = (
            (mode = :ka, ka_backend = KernelAbstractions.CPU()),
            (mode = :ka, ka_backend = JLBackend()),
        )
        for policy in policies
            @testset "simulation: merged=$merged, $policy" begin
                actual = simulate_reservoir(case; precond = :cprw, policy..., options...)
                @test length(actual.states) == length(result.states)
                for (state, ref) in zip(actual.states, result.states)
                    @test state[:Pressure] ≈ ref[:Pressure] rtol = 1e-7
                    @test state[:Saturations] ≈ ref[:Saturations] rtol = 1e-7
                end
            end
        end
    end
end
