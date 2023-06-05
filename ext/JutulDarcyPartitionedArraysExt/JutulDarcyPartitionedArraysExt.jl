module JutulDarcyPartitionedArraysExt
    using Jutul, JutulDarcy
    using PrecompileTools
    # Specific dependencies
    using PartitionedArrays, MPI, HYPRE
    using LinearAlgebra

    import Jutul: PArraySimulator, MPISimulator, PArrayExecutor
    import Jutul: DebugPArrayBackend, JuliaPArrayBackend, MPI_PArrayBackend
    import Jutul: partition_distributed, simulate_parray, @tic
    import JutulDarcy: reservoir_partition, partitioner_input

    function JutulDarcy.setup_reservoir_simulator_parray(
            case::JutulCase,
            backend::PArrayBackend;
            conn = :unit,
            np = missing,
            kwarg...
        )
        if ismissing(np)
            np = MPI.Comm_size(MPI.COMM_WORLD)
        end
        np::Int
        N, T, groups = partitioner_input(case.model, case.parameters, conn = conn)
        rmodel = reservoir_model(case.model)
        nc = number_of_cells(rmodel.domain)
        p_num = partition_distributed(N, T, nc = nc, np = np, groups = groups)
        p = reservoir_partition(case.model, p_num)
        return PArraySimulator(case, p; backend = backend, kwarg...)
    end

    function JutulDarcy.set_default_cnv_mb!(config::JutulConfig, sim::Jutul.PArraySimulator; kwarg...)
        simulators = sim.storage[:simulators]
        map(simulators, config[:configs]) do sim, cfg
            JutulDarcy.set_default_cnv_mb!(cfg, sim)
        end
        return config
    end

    function Jutul.parray_preconditioner_apply!(global_out, main_prec::CPRPreconditioner{<:BoomerAMGPreconditioner, <:Any}, X, preconditioners, simulator, arg...)
        global_cell_vector = simulator.storage.distributed_cell_buffer
        global_buf = simulator.storage.distributed_residual_buffer
        map(local_values(X), preconditioners, ghost_values(X)) do x, prec, x_g
            @. x_g = 0.0
            JutulDarcy.apply_cpr_first_stage!(prec, x, arg...)
            nothing
        end
        # The following is an unsafe version of this:
        # copy!(global_cell_vector, main_prec.p)
        p_h = main_prec.p
        @assert !isnothing(p_h) "CPR is not properly initialized."
        map(own_values(global_cell_vector), preconditioners) do ov, prec
            helper = prec.pressure_precond.data[:assembly_helper]
            indices = helper.indices
            indices::Vector{HYPRE.HYPRE_BigInt}
            nvalues = indices[end] - indices[1] + 1
            HYPRE.@check HYPRE.HYPRE_IJVectorGetValues(p_h, nvalues, indices, ov)
        end
        # End unsafe shenanigans

        # consistent!(global_cell_vector) |> wait
        map(own_values(global_buf), own_values(global_cell_vector), preconditioners) do dx, dp, prec
            bz = prec.block_size
            for i in eachindex(dp)
                JutulDarcy.set_dp!(dx, bz, dp, i)
            end
            nothing
        end
        if main_prec.full_system_correction
            mul_ix = nothing
        else
            mul_ix = 1
        end
        full_op = Jutul.parray_linear_system_operator(simulator.storage.simulators, global_buf)
        mul!(X, full_op, global_buf, -1.0, true)

        map(local_values(global_out), local_values(X), preconditioners, local_values(global_cell_vector), ghost_values(X)) do y, x, prec, dp, x_g
            @. x_g = 0.0
            apply!(y, prec.system_precond, x, arg...)
            bz = prec.block_size
            JutulDarcy.increment_pressure!(y, dp, bz)
            nothing
        end
        consistent!(global_out) |> wait
        global_out
    end
    @compile_workload begin
        targets = [(true, :csc), (true, :csr)]
        # MPI, trivial partition
        JutulDarcy.precompile_darcy_multimodels(targets,
            dims = (4, 1, 1),
            setuparg = (mode = :mpi, ),
            split_wells = true
        )
        # Native PArray, non-trivial partition
        JutulDarcy.precompile_darcy_multimodels(targets,
            dims = (4, 1, 1),
            default_linsolve = false,
            setuparg = (
                mode = :parray,
                parray_arg = (np = 2, ),
                precond = :ilu0
                ),
            split_wells = true
        )
    end
end
