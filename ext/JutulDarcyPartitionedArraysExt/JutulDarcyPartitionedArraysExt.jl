module JutulDarcyPartitionedArraysExt
    using Jutul, JutulDarcy
    # Specific dependencies
    using PartitionedArrays, MPI, HYPRE
    using LinearAlgebra

    import Jutul: PArraySimulator, MPISimulator, PArrayExecutor
    import Jutul: DebugPArrayBackend, JuliaPArrayBackend, MPI_PArrayBackend
    import Jutul: partition_distributed, simulate_parray
    import JutulDarcy: reservoir_partition, partitioner_input

    function JutulDarcy.setup_reservoir_simulator_parray(
            case::JutulCase,
            backend::PArrayBackend;
            conn = :unit,
            np = MPI.Comm_size(MPI.COMM_WORLD)
        )
        N, T, groups = partitioner_input(case.model, case.parameters, conn = conn)
        rmodel = reservoir_model(case.model)
        nc = number_of_cells(rmodel.domain)
        p_num = partition_distributed(N, T, nc = nc, np = np, groups = groups)
        p = reservoir_partition(case.model, p_num)
        return PArraySimulator(case, p, backend = backend)
    end

    function Jutul.parray_preconditioner_apply!(global_out, main_prec::CPRPreconditioner{<:BoomerAMGPreconditioner, <:Any}, X, preconditioners, simulator, arg...)
        tmr = simulator.storage.global_timer
        global_cell_vector = simulator.storage.distributed_cell_buffer
        global_buf = simulator.storage.distributed_residual_buffer
        # tic!(tmr)
        # toc!(tmr, "Apply preconditioner communication")
        # tic!(tmr)
        map(local_values(X), preconditioners, ghost_values(X)) do x, prec, x_g
            @. x_g = 0.0
            JutulDarcy.apply_cpr_first_stage!(prec, x, arg...)
            nothing
        end
        # toc!(tmr, "Apply CPR stage 1")
        # tic!(tmr)
        # The following is an unsafe version of this:
        # copy!(global_cell_vector, main_prec.p)
        map(own_values(global_cell_vector), preconditioners) do ov, prec
            helper = prec.pressure_precond.data[:assembly_helper]
            indices = helper.indices
            indices::Vector{HYPRE.HYPRE_BigInt}
            nvalues = indices[end] - indices[1] + 1
            HYPRE.@check HYPRE.HYPRE_IJVectorGetValues(main_prec.p, nvalues, indices, ov)
        end
        # End unsafe shenanigans

        # consistent!(global_cell_vector) |> wait
        # toc!(tmr, "Apply preconditioner communication")
        # tic!(tmr)

        map(own_values(global_buf), own_values(global_cell_vector), preconditioners) do dx, dp, prec
            bz = prec.block_size
            for i in eachindex(dp)
                JutulDarcy.set_dp!(dx, bz, dp, i)
            end
            nothing
        end
        # n = length(global_buf)
        if main_prec.full_system_correction
            mul_ix = nothing
        else
            mul_ix = 1
        end
        # distributed_mul! = Jutul.setup_parray_mul!(tmr, simulator.storage.simulators, mul_ix)
        full_op = Jutul.parray_linear_system_operator(tmr, simulator.storage.simulators, global_buf)
        # full_op = LinearOperator(Float64, n, n, false, false, distributed_mul!)
        mul!(X, full_op, global_buf, -1.0, true)

        map(local_values(global_out), local_values(X), preconditioners, local_values(global_cell_vector), ghost_values(X)) do y, x, prec, dp, x_g
            @. x_g = 0.0
            apply!(y, prec.system_precond, x, arg...)
            bz = prec.block_size
            JutulDarcy.increment_pressure!(y, dp, bz)
            nothing
        end
        # toc!(tmr, "CPR second stage")
        # tic!(tmr)
        consistent!(global_out) |> wait
        # toc!(tmr, "Apply preconditioner communication")
        global_out
    end

end
