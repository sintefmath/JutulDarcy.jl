function NLDDSimulator(model, partition = missing;
        state0 = setup_state(model),
        parameters = setup_parameters(model),
        kwarg...
    )
    case = JutulCase(model, parameters = parameters, state0 = state0)
    return NLDDSimulator(case, partition; kwarg...)
end

function NLDDSimulator(case::JutulCase, partition = missing;
        submodels = nothing,
        check_buffers = false,
        executor = Jutul.default_executor(),
        cells_per_block = missing,
        mpi_sync_after_solve = true,
        no_blocks = missing,
        kwarg...
    )
    (; model, state0, parameters) = case
    if ismissing(partition)
        rmodel = JutulDarcy.reservoir_model(model)
        G = rmodel.domain
        nc = number_of_cells(G)
        if ismissing(cells_per_block) && ismissing(no_blocks)
            cells_per_block = 1000
        end

        if ismissing(cells_per_block)
            @assert !ismissing(no_blocks)
            N = no_blocks
        else
            @assert ismissing(no_blocks)
            N = Int64(clamp(ceil(nc/cells_per_block), 1, nc))
        end
        p = partition_from_N(model, parameters, N)
        partition = reservoir_partition(model, p);
    elseif partition isa Vector{Int}
        # Convert it.
        partition = reservoir_partition(model, partition)
    end
    partition::Jutul.AbstractDomainPartition
    outer_sim = Simulator(model; state0 = state0, parameters = deepcopy(parameters), executor = executor, kwarg...)
    np = number_of_subdomains(partition)
    maximum_np = np
    if isnothing(submodels)
        is_distributed_solve = executor isa Jutul.PArrayExecutor
        if is_distributed_solve
            n_self = executor.data[:n_self]
            nc = length(executor.data[:partition])
            active_global = [i <= n_self for i in 1:nc]
            is_mpi = executor.mode isa Jutul.MPI_PArrayBackend
            if is_mpi
                maximum_np = Jutul.mpi_scalar_allreduce(np, max, executor)
            end
        else
            is_mpi = false
            active_global = missing
        end
        submodels = build_submodels(model, partition, active_global = active_global)

        if is_distributed_solve && false
            categorized = zeros(Int, nc)
            categorized_bnd = zeros(Int, nc)

            for m in submodels
                rm = reservoir_model(m)
                gmap = rm.domain.global_map
                gmap::Jutul.FiniteVolumeGlobalMap
                for (cell, is_bnd) in zip(gmap.cells, gmap.cell_is_boundary)
                    if is_bnd
                        categorized_bnd[cell] += 1
                    else
                        categorized[cell] += 1
                    end
                    if !active_global[cell]
                        @assert is_bnd
                    end
                end
            end
            for i in 1:nc
                if i > n_self
                    @assert categorized_bnd[i] > 0 "Cell $i was categorized $(categorized_bnd[i]) times as boundary"
                else
                    @assert categorized[i] == 1 "Cell $i was categorized as interior $(categorized[i]) times."
                end
            end
        end
    end
    state0 = setup_state(model, state0)
    subsim = Vector{Simulator}(undef, np)
    coarse_neighbors = Vector{Vector{Int}}()
    @showprogress "Building $np simulators " for i = 1:np
        sm = submodels[i]
        substate0 = substate(state0, model, sm, :variables)
        subparam = substate(parameters, model, sm, :parameters)
        subsim[i] = Simulator(sm, state0 = substate0, parameters = subparam)

        neigh_c = Jutul.coarse_neighborhood(partition, sm)
        filter!(!isequal(i), neigh_c)
        push!(coarse_neighbors, neigh_c)
    end
    storage = JutulStorage()
    buf = map(s -> primary_variable_buffer(s.storage, s.model), subsim)
    buf_pos = map(x -> primary_variable_positions(outer_sim, x[2], buf[x[1]], x[1], partition), enumerate(subsim))
    if check_buffers
        ok = check_primary_variable_positions(outer_sim.model, buf_pos)
        if !ok
            error("Something went wrong with primary variable positions")
        end
    end
    storage[:primary_variable_buffers] = buf
    storage[:primary_variable_positions] = buf_pos
    if reservoir_model(model) isa JutulDarcy.StandardBlackOilModel
        storage[:black_oil_primary_buffers] = map(s -> black_oil_primary_buffers(s.storage, s.model), subsim)
    end
    storage[:local_changes] = map(s -> reservoir_change_buffers(s.storage, s.model), subsim)
    storage[:boundary_discretizations] = map((s) -> boundary_discretization(outer_sim.storage, outer_sim.model, s.model, s.storage), subsim)
    storage[:solve_log] = NLDDSolveLog()
    storage[:coarse_neighbors] = coarse_neighbors
    storage[:maximum_number_of_subdomains_globally] = maximum_np
    storage[:is_mpi] = is_mpi
    storage[:mpi_sync_after_solve] = mpi_sync_after_solve

    # storage = convert_to_immutable_storage(storage)
    NLDDSimulator(outer_sim, executor, partition, nothing, subsim, storage)
end

@enum LocalSolveStatus begin
    local_solve_skipped
    local_already_converged
    local_solved_in_single_step
    local_solved_in_multiple_steps
    local_solve_failure
end

global_forces(forces) = forces.outer
local_forces(forces, i) = forces.subdomains[i]

function Jutul.perform_step!(
        simulator::NLDDSimulator, dt, forces, config;
        iteration = 0,
        solve = true,
        report = Jutul.setup_ministep_report(),
        executor = Jutul.default_executor(),
        kwarg...
    )
    log = simulator.storage[:solve_log]
    strategy = config[:strategy]

    base_arg = (simulator, dt, forces, config, iteration)
    if executor isa Jutul.DefaultExecutor
        # This means we are not in MPI and we must call this part manually. This
        # contains the entire local stage of the algorithm.
        Jutul.perform_step_per_process_initial_update!(
            simulator, dt, forces, config;
            iteration = iteration,
            report = report,
            executor = executor
        )
    end
    active = report[:local_solves_active]
    subreports = report[:subdomains]
    status = report[:solve_status]

    e = NaN
    converged = false
    @tic "global" try
        # Then perform step with the outer solver
        e, converged, report = global_stage(base_arg..., status, subreports;
            report = report,
            newton_fallback = !active,
            solve = solve,
            executor = executor,
            kwarg...
        )
        store_metadata_nldd_logger!(log, report, active)
    catch excptn
        if config[:failure_cuts_timestep]
            @warn "Exception occured in NLDD solve: $excptn. Attempting to cut time-step since failure_cuts_timestep = true."
            report[:failure_exception] = excptn
        else
            rethrow(excptn)
        end
    end
    return (e, converged, report)
end

function Jutul.initialize_before_first_timestep!(sim::NLDDSimulator, first_dT; forces = forces, config = config)
    @tic "solve" begin
        @tic "global" begin
            s = sim.simulator
            @tic "secondary variables" Jutul.update_secondary_variables!(s.storage, s.model)
        end
        @tic "local" begin
            @tic "secondary variables" for s in sim.subdomain_simulators
                Jutul.update_secondary_variables!(s.storage, s.model)
            end
        end
    end
end

function flag_as_done(solve_ok::Vector{Bool}, i)
    solve_ok[i] = true
end
function flag_as_done(solve_ok, i)
    solve_ok[i][] = true
end
function subdomain_ok(solve_ok, i)
    return solve_ok[i][]
end
function did_solve(subreports, i)
    r = subreports[i]
    solved = !isnothing(r) && (length(r) > 1 || length(r[1][:steps]) > 1)
    # if !solved
    #     @info "Subdomain $i skipped solving"
    # end
    return solved
end

function local_stage(simulator, dt, forces, config, iteration)
    is_aspen = config[:method] == :aspen
    is_gauss_seidel = config[:gauss_seidel] && !is_aspen
    # Perform DD pass (nonlinear solves)
    sim_global = simulator.simulator
    sub_sims = simulator.subdomain_simulators
    sim_order = missing
    n = length(sub_sims)

    configs = config[:config_subdomains]
    subreports = Vector{Any}(undef, n)
    fill!(subreports, nothing)

    use_threads = config[:nldd_threads]
    solve_status = fill(local_solve_skipped, n)

    @debug "Solving local systems..."
    failures = Vector{Int64}()

    do_solve = zeros(Bool, n)
    allow_early_termination = config[:subdomain_failure_cuts]

    function pre(i)
        subsim = sub_sims[i]
        change_buffer = simulator.storage.local_changes[i]
        store_reservoir_change_buffer!(change_buffer, subsim, config)
        # Make sure that recorder is updated to match global context
        rec_global = sim_global.storage.recorder
        rec_local = subsim.storage.recorder
        Jutul.reset!(rec_local, rec_global)
        # Update state0
        # state0 should have the right secondary variables since it comes from a converged state.
        update_subdomain_from_global(sim_global, subsim, i, current = false, transfer_secondary = true)
        # We do not transfer the secondaries. This means that they are
        # updated/recomputed locally. The current primary variables in the
        # global state have been updated by the last global stage and the
        # secondary variables are then "out of sync". The local solvers do the
        # work of recomputing them.
        update_subdomain_from_global(sim_global, subsim, i, current = true, transfer_secondary = iteration == 1)
        should_solve = check_if_subdomain_needs_solving(change_buffer, subsim, config, iteration)
        do_solve[i] = should_solve
    end
    function solve_and_transfer(i)
        should_solve = do_solve[i]
        sim = sub_sims[i]
        cfg = configs[i]
        il_i = cfg[:info_level]
        if should_solve
            if il_i > 0
                jutul_message("Subdomain: $i", "Solving subdomain.")
            end
            # cfg[:always_update_secondary] = true
            ok, rep = solve_subdomain(sim, i, dt, forces, cfg)
        else
            if il_i > 0
                jutul_message("Subdomain: $i", "Skipping subdomain.")
            end
            ok = true
            rep = nothing
        end
        subreports[i] = rep
        solve_status[i] = get_solve_status(rep, ok)
        # Still need to transfer over secondary variables
        if !ok
            if config[:info_level] > 0
                Jutul.jutul_message("Failure $i", color = :red)
            end
            update_subdomain_from_global(sim_global, sim, i, current = true, transfer_secondary = false)
        end
        update_global_from_subdomain(sim_global, sim, i, secondary = true)
        return ok
    end
    # Now that we have the main functions set up we can do either kind of local solves
    t = @elapsed if is_gauss_seidel
        # Gauss-seidel (non-deterministic if threads are on)
        function gs_solve(i)
            pre(i)
            return solve_and_transfer(i)
        end
        @tic "local solves [GS]" sim_order = gauss_seidel_for_each_subdomain_do(
            gs_solve,
            simulator,
            sub_sims,
            subreports,
            config[:gauss_seidel_order],
            allow_early_termination
        )
    else
        execute = (f) -> for_each_subdomain_do(f, n, use_threads, config[:nldd_thread_type], allow_early_termination)
        # Jacobi solve
        @tic "prepare local solves" execute(pre)
        if is_aspen
            # Need to store primary variables before any local solves.
            @tic "store primary" store_subdomain_primary_variables!(simulator, config, do_delta = false)
        end
        @tic "local solves [Jacobi]" execute(solve_and_transfer)
    end
    for (i, st) in enumerate(solve_status)
        if st == local_solve_failure
            push!(failures, i)
        end
    end
    bad = length(failures)
    if bad > 0
        msg = "It $iteration: $bad of $n subdomains failed to converge"
        if config[:subdomain_failure_throws]
            error(msg)
        elseif config[:info_level] > 0
            if config[:info_level] == 1
                Jutul.jutul_message("Convergence", msg, color = :light_yellow)
            else
                @warn "$msg:" failures
                for i in eachindex(subreports)
                    if solve_status[i] == local_solve_failure
                        sr = subreports[i][end][:steps]
                        f_r = sr[end]
                        nit = length(sr)
                        if haskey(f_r, :linear_iterations)
                            l_its = f_r[:linear_iterations]
                        else
                            l_its = NaN
                        end
                        Jutul.jutul_message("Ω$i", "$(l_its) linear iterations, $nit iterations:")
                        Jutul.get_convergence_table(f_r[:errors], 3, nit, config)
                    end
                end
            end
        end
    end
    if config[:info_level] > 1
        m = sum(x-> !(x == local_solve_skipped), solve_status)
        if m > 0
            tot_str = Jutul.get_tstr(t, 1)
            avg_str = Jutul.get_tstr(t/m, 1)
            jutul_message("NLDD", "Solved $m/$n domains in $tot_str ($avg_str average)")
        else
            jutul_message("NLDD", "Solved no subdomains.")
        end
    end
    return subreports, sim_order, t, solve_status, failures
end

function get_solve_status(rep, ok)
    if isnothing(rep)
        status = local_solve_skipped
    else
        if ok
            if length(rep) > 1
                status = local_solved_in_multiple_steps
            else
                rep = only(rep)
                nstep = length(rep[:steps])
                if nstep == 0
                    status = local_already_converged
                elseif nstep == 1 && !haskey(rep[:steps][1], :linear_solver)
                    status = local_already_converged
                else
                    status = local_solved_in_single_step
                end
            end
        else
            status = local_solve_failure
        end
    end
    return status
end

function for_each_subdomain_do(f, n, use_threads, thread_type, early_stop::Bool)
    if use_threads
        if thread_type == :batch
            @batch for i = 1:n
                ok_i = f(i)
                if early_stop && !ok_i
                    break
                end
            end
        elseif thread_type == :default
            disable_polyester_threads() do
                Threads.@threads for i = 1:n
                    ok_i = f(i)
                    if early_stop && !ok_i
                        break
                    end
                end
            end
        else
            disable_polyester_threads() do
                if thread_type == :fair
                    worksteal_par!(f, n, Threads.FairSchedule())
                else
                    worksteal_par!(f, n, Threads.StaticSchedule())
                end
            end
        end
    else
        for i = 1:n
            ok_i = f(i)
            if early_stop && !ok_i
                break
            end
        end
    end
end

function gauss_seidel_for_each_subdomain_do(f, sim, simulators, subreports, strategy::Symbol, early_stop::Bool)
    n = length(simulators)
    is_mpi = sim.storage.is_mpi
    should_sync = sim.storage.mpi_sync_after_solve
    has_sync = haskey(sim.executor.data, :distributed_primary_variables_sync_function)
    if is_mpi && has_sync && should_sync
        sync_function = sim.executor.data[:distributed_primary_variables_sync_function]
        n_total = sim.storage[:maximum_number_of_subdomains_globally]
    else
        sync_function = (; kwarg...) -> nothing
        n_total = n
    end
    num_solved = 0
    function solve_gauss_seidel_iteration!(i)
        f(i)
        sync_function()
        num_solved += 1
    end
    if strategy == :adaptive
        sim_order = Int[]
        coarse_neighbors = sim.storage.coarse_neighbors
        sizehint!(sim_order, n)
        injectors = find_injectors(simulators)
        if length(injectors) == 0
            v, highest = findmax(i -> find_max_interior_pressure(simulators[i]), eachindex(simulators))
            injectors = [highest]
        end
        num_neighbor_updates = zeros(Int, n)
        n_solved = 0

        function solve_adaptive(i)
            solve_gauss_seidel_iteration!(i)
            push!(sim_order, i)
            subreports_i = subreports[i]
            if !isnothing(subreports_i)
                n_its = 0
                for rep in subreports_i
                    n_its += length(rep[:steps])
                end
                for j in coarse_neighbors[i]
                    if num_neighbor_updates[j] != -1
                        num_neighbor_updates[j] += n_its
                    end
                end
            end
            num_neighbor_updates[i] = -1
            n_solved += 1
        end
        # Solve injectors first
        for i in injectors
            solve_adaptive(i)
        end
        # Then solve domains in order of adjacent local solves
        while n_solved < n
            v, ix = findmax(num_neighbor_updates)
            @assert v != -1
            solve_adaptive(ix)
        end
        # @assert length(unique(sim_order)) == n
    else
        if strategy == :linear
            sim_order = eachindex(simulators)
        else
            if strategy == :pressure
                sort_function = x -> -find_max_interior_pressure(x)
            elseif strategy == :potential
                sort_function = x -> -find_block_potential(x)
            else
                error("Ordering $strategy is not supported.")
            end
            pval = map(sort_function, simulators)
            sim_order = sortperm(pval)
        end
        for i in sim_order
            ok_i = solve_gauss_seidel_iteration!(i)
            if early_stop && !ok_i
                break
            end
        end
    end
    if is_mpi
        for i in (num_solved+1):n_total
            # Need to account for other process waiting on a (potentially trivial) sync.
            sync_function()
        end
    end
    return sim_order
end

function find_injectors(simulators)
    injectors = Int[]
    for i in eachindex(simulators)
        sim = simulators[i]
        model = sim.model
        if model isa Jutul.SimulationModel
            continue
        end
        for (k, v) in pairs(model.models)
            if v.domain isa JutulDarcy.WellControllerDomain
                cfg = sim.storage.state[k].WellGroupConfiguration
                for (wk, wv) in pairs(cfg.operating_controls)
                    if wv isa JutulDarcy.InjectorControl
                        push!(injectors, i)
                    end
                end
            end
        end
    end
    # Sort so that highest pressure comes first
    pval = map(find_max_interior_pressure, simulators)
    sort!(injectors, by = i -> -pval[i])
    # sort!(injectors, by = i -> -find_max_interior_pressure(simulators[i]))
    return unique!(injectors)
end

function find_max_interior_pressure(sim)
    return find_max_interior_pressure(sim.model, sim.storage.state)
end

function find_max_interior_pressure(model::MultiModel, state)
    p_max = -Inf
    for k in Jutul.submodel_symbols(model)
        model_k = model[k]
        state_k = state[k]
        if haskey(state_k, :Pressure)#  && k == :Reservoir
            mm = global_map(model_k)::Union{Jutul.FiniteVolumeGlobalMap{Int}, Jutul.TrivialGlobalMap}
            p_max_k = find_max_interior_pressure(model_k, state_k, mm)::Float64
            p_max = max(p_max, p_max_k)::Float64
        end
    end
    return p_max
end

function find_max_interior_pressure(model, state)
    return find_max_interior_pressure(model, state, global_map(model.domain))
end

function find_max_interior_pressure(model, state, g_map::Jutul.FiniteVolumeGlobalMap)
    is_bnd = g_map.cell_is_boundary
    p = state.Pressure
    return find_max_pressure(is_bnd, p)
end

function find_max_interior_pressure(model, state, ::Jutul.TrivialGlobalMap)
    p = state.Pressure
    return maximum(value, p)::Float64
end

function find_max_pressure(is_bnd, p)
    p_max = -Inf
    for i in eachindex(p)
        if !is_bnd[i]
            p_max = max(p_max, value(p[i]))::Float64
        end
    end
    return p_max::Float64
end

function find_block_potential(sim)
    rmodel = reservoir_model(sim.model)
    rstate = sim.storage.Reservoir.state
    is_bnd = rmodel.domain.global_map.cell_is_boundary
    cc = rmodel.data_domain[:cell_centroids]
    p = rstate.Pressure
    if size(cc, 1) == 3
        z = view(cc, 3, :)
        pv = rstate.FluidVolume
        rho = rstate.PhaseMassDensities
        S = rstate.Saturations
        find_block_potential(is_bnd, p, rho, S, pv, z)
    else
        return find_max_pressure(is_bnd, p)
    end
end

function find_block_potential(is_bnd, p, rho, S, pv, z, g = Jutul.gravity_constant)
    rho_avg = 0.0
    pv_t = 0.0
    z_avg = 0.0
    p_avg = 0.0
    for i in eachindex(is_bnd)
        if !is_bnd[i]
            pv_i = value(pv[i])
            rho_i = 0.0
            for j in axes(S, 1)
                rho_i += value(S[j, i])*value(rho[j, i])
            end
            p_avg += value(p[i])*pv_i
            rho_avg += rho_i*pv_i
            z_avg += value(z[i])*pv_i
            pv_t += pv_i
        end
    end
    z_avg /= pv_t
    rho_avg /= pv_t
    p_avg /= pv_t
    # Finally compute potential
    pot = p_avg + rho_avg*z_avg*g
    return pot
end

function worksteal_par!(solve, n, s = Threads.StaticSchedule())
    chnl = Channel() do ch
        foreach(i -> put!(ch, i), 1:n)
    end
    Threads.foreach(solve, chnl, schedule=s)
end

function global_stage(
        simulator, dt, forces, config, iteration, solve_status, subreports;
        report = Jutul.setup_ministep_report(),
        newton_fallback = false,
        solve = true,
        sim_kwarg...
    )
    # @warn "Starting global stage"
    m = config[:method]
    g_forces = global_forces(forces)
    s = simulator.simulator
    function solve_fi(do_solve = solve)
        # Secondary variables get copied over - no need to update them.
        return Jutul.perform_step!(
            s, dt, g_forces, config;
            report = report,
            iteration = iteration,
            update_secondary = newton_fallback,
            solve = do_solve,
            sim_kwarg...
        )
    end
    # If all subdomains converged, we can return early
    all_local_converged = true
    any_failures = false
    if isnothing(subreports)
        all_local_converged = false
    else
        for i in eachindex(subreports)
            status_i = solve_status[i]
            ok_i = status_i == local_already_converged
            any_failures = any_failures || status_i == local_solve_failure
            all_local_converged = all_local_converged && ok_i
        end
    end
    if all_local_converged && (iteration > config[:min_nonlinear_iterations]) && config[:subdomain_tol_sufficient]
        report[:secondary_time] = 0.0
        report[:equations_time] = 0.0
        report[:linear_system_time] = 0.0
        report[:converged] = true
        report[:convergence_time] = 0.0
        Jutul.extra_debug_output!(report, s.storage, s.model, config, iteration, dt)
        out = (0.0, true, report)
        il = config[:info_level]
        if il > 1
            errors = Dict()
            Jutul.get_convergence_table(errors, il, iteration, config)
        end
    elseif config[:subdomain_failure_cuts] && any_failures
        report[:failure] = true
        report[:failure_exception] = "Subdomains failed and subdomain_failure_cuts = true"
        return (1.0, false, report)
    else
        if m == :nldd
            n = length(simulator.subdomain_simulators)
            for i in 1:n
                sub_sim = simulator.subdomain_simulators[i]
                if solve_status[i] == local_solve_failure
                    # Jutul.jutul_message("Failure $i", color = :red)
                    # reset_local_domain!(simulator.storage, s, sub_sim, forces, dt, i)
                end
            end
            out = solve_fi()
        elseif m == :aspen
            if newton_fallback
                out = solve_fi()
            else
                out = global_aspen!(simulator, dt, forces, g_forces, config, iteration, solve_status, subreports; report = report, sim_kwarg...)
            end
        elseif m == :local_only
            out = solve_fi(false)
        else
            error("Method must be either :aspen or :nldd. Method was: :$m")
        end
    end
    check_locals_after = config[:debug_checks]
    if check_locals_after && out[2] == true
        if simulator.simulator.storage.LinearizedSystem isa MultiLinearizedSystem
            @info "Outer converged, checking inner..."
            sub_sims = simulator.subdomain_simulators
            n = length(sub_sims)
            configs = config[:config_subdomains]
            for i in eachindex(sub_sims)
                sim, cfg = sub_sims[i], configs[i]
                cfg[:info_level] = 3
                # Update state0
                update_subdomain_from_global(simulator.simulator, sub_sims[i], i, current = false)
                # Update state
                update_subdomain_from_global(simulator.simulator, sub_sims[i], i, current = true, transfer_secondary = false)
                ok, subrep = solve_subdomain(sim, i, dt, forces, cfg)
                @assert ok
                nsubstep = length(subrep[1][:steps]) - 1
                @assert nsubstep == 0 "Expected no substeps, was $nsubstep"
                # @info nsubstep cfg[:min_nonlinear_iterations]
                cfg[:info_level] = -1
            end
        end
    end
    return out
end


function reset_local_domain!(s, outer_sim, sub_sim, forces, dt, i)
    update_subdomain_from_global(outer_sim, sub_sim, i, current = true, transfer_secondary = true)
    s_i = sub_sim.storage
    m_i = sub_sim.model
    f_i = local_forces(forces, i)

    update_subdomain_from_global(outer_sim, sub_sim, i, current = false, transfer_secondary = true)
end


function solve_subdomain(sim, i, dt, forces, cfg)
    max_iter = cfg[:max_nonlinear_iterations]::Integer
    f = local_forces(forces, i)
    report = nothing
    ok = false

    try
        ok, report = Jutul.solve_timestep!(sim, dt, f, max_iter, cfg,
            finalize = true, prepare = true, update_explicit = false)
    catch e
        if e isa InterruptException
            # Don't leave the user trapped...
            rethrow(e)
        end
        @warn "Subdomain $i threw exception: $e"
    end
    return (ok, report)
end

function Jutul.update_before_step!(sim::NLDDSimulator, dt, forces; kwarg...)
    Jutul.update_before_step!(sim.simulator, dt, global_forces(forces); kwarg...)
    for (i, sim_i) in enumerate(sim.subdomain_simulators)
        Jutul.update_before_step!(sim_i, dt, local_forces(forces, i); kwarg...)
    end
end

function Jutul.update_after_step!(sim::NLDDSimulator, dt, forces; kwarg...)
    Jutul.update_after_step!(sim.simulator, dt, global_forces(forces); kwarg...)
end

function Jutul.reset_state_to_previous_state!(sim::NLDDSimulator)
    Jutul.reset_state_to_previous_state!(sim.simulator)
end

function Jutul.reset_previous_state!(sim::NLDDSimulator, state0)
    # It is sufficient to reset the "outer" simulator since everything
    # is copied over to local subdomains from that state.
    Jutul.reset_previous_state!(sim.simulator, state0)
end

function Jutul.store_output!(states, reports, step, sim::NLDDSimulator, config, report; kwarg...)
    return Jutul.store_output!(states, reports, step, sim.simulator, config, report; kwarg...)
end

function Jutul.select_nonlinear_relaxation(sim::NLDDSimulator, rel_type, reports, relaxation)
    return Jutul.select_nonlinear_relaxation(sim.simulator, rel_type, reports, relaxation)
end

function check_primary_variable_positions(model::MultiModel, positions)
    ok = true
    groups = model.groups
    if isnothing(groups)
        groups = ones(Int, length(model.models))
    end
    for current_group in unique(groups)
        num = 0
        for (group, m) in zip(groups, model.models)
            if group == current_group
                num += Jutul.number_of_degrees_of_freedom(m)
            end
        end
        count = zeros(Int, num)
        for pos in positions
            if length(pos) >= current_group
                for i in pos[current_group]
                    count[i] += 1
                end
            end
        end
        ok = ok && all(count .== 1)
    end
    return ok
end


function debug_comparision(sim_g, sim, i)
    @info "Checking $i" keys(sim.model.models)
    m = JutulDarcy.reservoir_model(sim.model)
    m_g = JutulDarcy.reservoir_model(sim_g.model)
    cells = m.domain.global_map.cells
    state = sim.storage.Reservoir.state
    state_g = sim_g.storage.Reservoir.state

    p = state[:Pressure]
    p_g = state_g[:Pressure][cells]
    @info "P" p p_g norm(p-p_g)/norm(p)

    S = state[:Saturations]
    S_g = state_g[:Saturations][:, cells]
    @info "P" S S_g norm(S-S_g)/norm(S)

    @warn "Before"
    for k in [:PhaseMassDensities, :TotalMasses, :RelativePermeabilities, :PhaseMobilities, :PhaseMassMobilities]
        x = state[k]
        x_g = state_g[k][:, cells]
        @info "$k" x x_g norm(x-x_g)/norm(x)
    end
    Jutul.update_secondary_variables!(sim_g.storage, sim_g.model)
    @warn "After"
    for k in [:PhaseMassDensities, :TotalMasses, :RelativePermeabilities, :PhaseMobilities, :PhaseMassMobilities]
        x = state[k]
        x_g = state_g[k][:, cells]
        @info "$k" x x_g norm(x-x_g)/norm(x)
    end
    kr = m[:RelativePermeabilities]
    kr_g = m_g[:RelativePermeabilities]
    @info "?" typeof(kr) typeof(kr_g) typeof(kr.regions) typeof(kr_g.regions)
    r = kr.regions
    r_g = kr_g.regions[cells]
    @error "Region check" r r_g norm(r - r_g) unique(r)
    for i in eachindex(r)
        reg = JutulDarcy.region(r, i)
        @assert reg == r_g[i]
    end

    c = length(cells)
    c_g = cells[end]
    kr_val = state[:RelativePermeabilities]
    kr_val_g = state_g[:RelativePermeabilities]
    indices = JutulDarcy.phase_indices(m.system)
    @info "Before" kr_val_g[:, c_g] kr_val[:, c]

    JutulDarcy.two_phase_relperm!(kr_val, S, kr.regions, kr.krw, kr.krow, indices, c)
    JutulDarcy.two_phase_relperm!(kr_val_g, state_g[:Saturations], kr_g.regions, kr_g.krw, kr_g.krow, indices, c_g)

    @info "Huh" kr_val_g[:, c_g] kr_val[:, c] length(kr_val)
    Jutul.update_secondary_variables!(sim.storage, sim.model)
    kr_val = sim.storage.Reservoir.state[:RelativePermeabilities]

    @info "Final" kr_val_g[:, c_g] kr_val[:, c] length(kr_val)


    error()
end

function Jutul.perform_step_per_process_initial_update!(simulator::NLDDSimulator, dt, forces, config;
        executor = Jutul.default_executor(),
        update_secondary = nothing,
        iteration = 0,
        report = Jutul.setup_ministep_report()
    )
    log = simulator.storage[:solve_log]
    strategy = config[:strategy]
    e = NaN
    converged = false

    base_arg = (simulator, dt, forces, config, iteration)
    if iteration == 1
        active = active_status_first_iteration(strategy, log)
        reset_nldd_logger!(log, active, forces, dt)
    else
        active = active_status(strategy, log, iteration)
    end
    @tic "local" if (active && config[:solve_subdomains])
        subreports, solve_order, t_sub, status, failures = local_stage(base_arg...)
    else
        subreports = nothing
        solve_order = nothing
        t_sub = 0.0
        failures = []
        status = [local_solve_skipped for i in eachindex(simulator.subdomain_simulators)]
    end
    report[:subdomains] = subreports
    report[:time_subdomains] = t_sub
    report[:subdomain_failures] = failures
    report[:local_solves_active] = active
    report[:solve_order] = solve_order
    report[:solve_status] = status

    return report
end
