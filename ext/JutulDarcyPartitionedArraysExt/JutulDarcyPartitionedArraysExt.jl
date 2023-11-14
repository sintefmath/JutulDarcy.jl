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

    function Jutul.parray_preconditioner_apply!(global_out, main_prec::CPRPreconditioner{<:BoomerAMGPreconditioner, <:Any}, R, preconditioners, simulator, arg...)
        global_cell_vector = simulator.storage.distributed_cell_buffer
        global_buf = simulator.storage.distributed_residual_buffer
        A_ps = main_prec.storage.A_ps

        @. global_out = 0.0
        npre = main_prec.npre
        npost = main_prec.npost
        if npre > 0
            JutulDarcy.apply_cpr_smoother!(global_out, R, global_buf, preconditioners, A_ps, npre)
        end
        @tic "cpr first stage" map(local_values(R), preconditioners, ghost_values(R)) do r, prec, x_g
            @. x_g = 0.0
            JutulDarcy.apply_cpr_pressure_stage!(prec, prec.storage, r, arg...)
            nothing
        end
        # The following is an unsafe version of this:
        # copy!(global_cell_vector, main_prec.p)
        p_h = main_prec.storage.p
        @assert !isnothing(p_h) "CPR is not properly initialized."
        @tic "hypre GetValues" map(own_values(global_cell_vector), preconditioners, own_values(global_out), own_values(global_buf)) do ov, prec, x_final, buf
            helper = prec.pressure_precond.data[:assembly_helper]
            bz = prec.storage.block_size
            indices = helper.indices
            indices::Vector{HYPRE.HYPRE_BigInt}
            nvalues = indices[end] - indices[1] + 1
            HYPRE.@check HYPRE.HYPRE_IJVectorGetValues(p_h, nvalues, indices, ov)
            JutulDarcy.increment_pressure!(x_final, ov, bz)
            if npost > 0
                @. buf = 0
                for i in 1:length(ov)
                    p_i = ov[i]
                    buf[(i-1)*bz + 1] = p_i
                end
            end
        end
        # End unsafe shenanigans
        if npost > 0
            JutulDarcy.correct_residual!(global_out, A_ps, global_buf)
            JutulDarcy.apply_cpr_smoother!(global_out, R, global_buf, preconditioners, A_ps, npost)
        end
        @tic "communication" consistent!(global_out) |> wait
        return global_out
    end

    function JutulDarcy.apply_cpr_smoother!(X::PVector, R::PVector, Buf::PVector, prec, A_ps, n; skip_last = false)
        for i in 1:n
            map(local_values(Buf), local_values(R), local_values(X), prec) do buf, r, x, p
                apply!(buf, p.system_precond, r)
                @. x += buf
            end
            if i < n || !skip_last
                JutulDarcy.correct_residual!(R, A_ps, X)
            end
        end
    end

    function Jutul.parray_update_preconditioners!(sim::Jutul.PArraySimulator, cpr::CPRPreconditioner{<:BoomerAMGPreconditioner, <:Any}, preconditioners, recorder)
        offset = sim.storage.process_offset
        n = sim.storage.nc_process
        comm = sim.storage.comm
        if sim.storage[:number_of_processes] > 1
            @assert sim.backend isa Jutul.MPI_PArrayBackend "Cannot use HYPRE with emulated multiple processes. Backend was $(sim.backend)"
        end

        function create_hypre_vector()
            x = HYPREVector(comm, offset + 1, offset + n)
            asm = HYPRE.start_assemble!(x)
            HYPRE.finish_assemble!(asm)
            return x
        end
        if isnothing(cpr.storage)
            A_p = HYPREMatrix(comm, offset + 1, offset + n)
            r_p = create_hypre_vector()
            p = create_hypre_vector()

            global_sol_buf = sim.storage.distributed_solution_buffer
            global_res_buf = sim.storage.distributed_residual_buffer
            A_ps = Jutul.parray_linear_system_operator(sim.storage.simulators, length(global_res_buf))
            p_sys = (A_p, r_p, p)
            rmodel = reservoir_model(sim.storage.model)
            bz = degrees_of_freedom_per_entity(rmodel, Cells())
            cpr.storage = JutulDarcy.CPRStorage(n, bz, A_ps, p_sys, global_sol_buf, global_res_buf)
        end
        cpr_storage = cpr.storage
        A_p = cpr_storage.A_p
        A_ps = cpr_storage.A_ps
        r_p = cpr_storage.r_p
        x_p = cpr_storage.p
        bz = cpr_storage.block_size
        w_rhs = cpr_storage.w_rhs

        map(sim.storage.simulators, preconditioners) do sim, prec
            storage = Jutul.get_simulator_storage(sim)
            model = Jutul.get_simulator_model(sim)
            sys = storage.LinearizedSystem
            prec.pressure_precond.data[:hypre_system] = (A_p, r_p, x_p)

            nc = number_of_cells(reservoir_model(model).domain)
            w_p = zeros(bz, nc)
            prec.storage = JutulDarcy.CPRStorage(A_p, r_p, x_p, missing, missing, A_ps, w_p, w_rhs, n, bz)
            Jutul.update_preconditioner!(prec, sys, model, storage, recorder, sim.executor)
            prec
        end
        return (cpr, preconditioners)
    end

    # @compile_workload begin
    #     targets = [(true, :csc), (true, :csr)]
    #     # MPI, trivial partition
    #     JutulDarcy.precompile_darcy_multimodels(targets,
    #         dims = (4, 1, 1),
    #         default_linsolve = false,
    #         setuparg = (
    #             mode = :mpi,
    #             precond = :ilu0
    #             ),
    #         split_wells = true
    #     )
    #     # Native PArray, non-trivial partition
    #     JutulDarcy.precompile_darcy_multimodels(targets,
    #         dims = (4, 1, 1),
    #         default_linsolve = false,
    #         setuparg = (
    #             mode = :parray,
    #             parray_arg = (np = 2, ),
    #             precond = :ilu0
    #             ),
    #         split_wells = true
    #     )
    # end
end
