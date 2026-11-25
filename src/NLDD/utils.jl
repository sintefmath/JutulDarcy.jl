import Krylov
export reservoir_partition, build_submodels

import JutulDarcy: set_default_cnv_mb!

function simulator_config(sim::NLDDSimulator;
        method = :nldd,
        subdomain_info_level = -1,
        inner_tol_mul = 1.0,
        inner_tol_final = 10.0,
        inner_max_timestep_cuts = 2,
        inner_min_nonlinear_iterations = 0,
        inner_max_nonlinear_iterations = missing,
        inner_accept_unconverged = false,
        inner_method = :nldd,
        tol_cnv = 1e-3,
        tol_mb = 1e-7,
        tol_cnv_well = 1e-2,
        tol_mb_well = 1e-3,
        tol_dp_well = 1e-3,
        inc_tol_dp_abs = Inf,
        inc_tol_dp_rel = Inf,
        inc_tol_dz = Inf,
        use_julia_amg = false,
        subdomain_precond = :ilu0,
        kwarg...
    )
    if subdomain_info_level isa Real
        slvl = subdomain_info_level
        subdomain_info_level = i -> slvl
    end
    check_before_solve = isinf(inc_tol_dp_abs) && isinf(inc_tol_dp_rel) && isinf(inc_tol_dz)
    cfg = simulator_config(sim.simulator)
    set_default_cnv_mb!(cfg, sim.simulator.model,
        tol_cnv = tol_cnv,
        tol_mb = tol_mb,
        tol_cnv_well = tol_cnv_well,
        tol_mb_well = tol_mb_well,
        tol_dp_well = tol_dp_well,
        inc_tol_dp_abs = inc_tol_dp_abs,
        inc_tol_dp_rel = inc_tol_dp_rel,
        inc_tol_dz = inc_tol_dz
    )
    is_mpi = sim.storage.is_mpi

    # Extra options
    add_option!(cfg, :nldd_threads, false, "Use threads for local solves. Not compatible with Gauss-Seidel.", types = Bool)
    add_option!(cfg, :nldd_thread_type, :default, "Type of threads to use", types = Symbol)
    add_option!(cfg, :gauss_seidel, method == :nldd, "Use Gauss-Seidel/multiplicative for local solves. Not compatible with nldd_threads = true", types = Bool)
    add_option!(cfg,
        :gauss_seidel_order,
        :pressure,
        "Order for Gauss-Seidel domain updates",
        types = Symbol,
        values = [:linear, :pressure, :potential, :adaptive]
    )
    add_option!(cfg, :nldd_max_sweeps, 1, "Maximum number of NLDD sweeps", types = Int)
    add_option!(cfg, :nldd_subdomain_groups, nothing, "Group of each subdomain. The local stage is then Gauss-Seidel within each group and Jacobi between groups (as with MPI, where each group is a process).", types = Union{Nothing, Vector{Int}})
    add_option!(cfg, :nldd_initial_newton_iterations, 0, "Number of global Newton iterations at the start of each time step before local solves are used.", types = Int)
    add_option!(cfg, :method, method, "Method to use", values = [:nldd, :aspen, :local_only], types = Symbol)
    add_option!(cfg, :debug_checks, false, "Enable extra expensive checks", types = Bool)
    add_option!(cfg, :solve_subdomains, true, "Solve subdomains", types = Bool)
    add_option!(cfg, :aspen_full_increment, false, "Solve full ASPEN update", types = Bool)
    add_option!(cfg, :aspen_increment_limits, true, "Apply the maximum relative/absolute changes of primary variables to the global ASPEN update. Bounds on the values (e.g. saturations in [0, 1]) are always applied.", types = Bool)
    add_option!(cfg, :strategy, DefaultNLDDStrategy(), "Strategy to use for applying NLDD/ASPEN")
    add_option!(cfg, :inner_tol_mul, inner_tol_mul, "Factor applied to the convergence tolerances of the subdomains relative to the global tolerances. Either a number (applied to all tolerances) or a NamedTuple with tolerance names as keys (e.g. (tol_mb = 100.0, tol_cnv = 0.1)), where tolerances that are not given keep the global value.", types = Union{Real, NamedTuple})
    # Only subdomains that converge without iterating count towards global
    # convergence, so inner_tol_final (only used at the last iteration) does not
    # matter here.
    same_tol = all(x -> x <= 1.0, inner_tol_mul isa Real ? (inner_tol_mul,) : values(inner_tol_mul))
    add_option!(cfg, :subdomain_tol_sufficient, same_tol && !is_mpi, "Tolerances in subdomains are at least tight enough to be able to conclude global convergence.", types = Bool)

    # Subdomain tolerances for when to solve a local subdomain
    add_option!(cfg, :solve_tol_temperature, nothing, "Local subdomains are solved if maximum temperature change at boundary exceeds this value.", types = Union{Float64, Nothing})
    add_option!(cfg, :solve_tol_temperature_mean, nothing, "Local subdomains are solved if mean of temperature change at boundary exceeds this value.", types = Union{Float64, Nothing})

    add_option!(cfg, :solve_tol_pressure, nothing, "Local subdomains are solved if maximum pressure change at boundary exceeds this value.", types = Union{Float64, Nothing})
    add_option!(cfg, :solve_tol_pressure_mean, nothing, "Local subdomains are solved if mean of pressure change at boundary exceeds this value.", types = Union{Float64, Nothing})

    add_option!(cfg, :solve_tol_saturations, nothing, "Local subdomains are solved if maximum saturation change at boundary exceeds this value.", types = Union{Float64, Nothing})
    add_option!(cfg, :solve_tol_saturations_mean, nothing, "Local subdomains are solved if mean of saturation change at boundary exceeds this value.", types = Union{Float64, Nothing})

    add_option!(cfg, :solve_tol_densities, nothing, "Local subdomains are solved if maximum density change at boundary exceeds this value.", types = Union{Float64, Nothing})
    add_option!(cfg, :solve_tol_densities_mean, nothing, "Local subdomains are solved if mean of density change at boundary exceeds this value.", types = Union{Float64, Nothing})

    add_option!(cfg, :solve_tol_mobility, nothing, "Local subdomains are solved if maximum mobility at boundary (relative to to total mobility) exceeds this value.", types = Union{Float64, Nothing})
    add_option!(cfg, :solve_tol_mobility_mean, nothing, "Local subdomains are solved if mean of mobility at boundary (relative to to total mobility) exceeds this value.", types = Union{Float64, Nothing})

    add_option!(cfg, :solve_tol_composition, nothing, "Local subdomains are solved if maximum change in composition at boundary exceeds this value.", types = Union{Float64, Nothing})
    add_option!(cfg, :solve_tol_composition_mean, nothing, "Local subdomains are solved if mean of change in composition at boundary exceeds this value.", types = Union{Float64, Nothing})

    add_option!(cfg, :solve_tol_phase_mass_fractions, nothing, "Local subdomains are solved if maximum change in phase mass fractions at boundary exceeds this value.", types = Union{Float64, Nothing})
    add_option!(cfg, :solve_tol_phase_mass_fractions_mean, nothing, "Local subdomains are solved if mean of change in phase mass fractions at boundary exceeds this value.", types = Union{Float64, Nothing})

    add_option!(cfg, :solve_tol_first_newton, true, "Skip local solves for first iteration when tolerances are used.", types = Bool)

    add_option!(cfg, :subdomain_failure_throws, false, "Throw exception upon subdomain failure to solve.", types = Bool)
    add_option!(cfg, :subdomain_failure_cuts, false, "Cut outer timestep if any subdomain fails to solve.", types = Bool)

    add_option!(cfg, :always_solve_wells, false, "Always solve wells in local subdomain", types = Bool)

    # Subdomain setup
    subsims = sim.subdomain_simulators
    n = length(subsims)
    subconfigs = Vector{JutulConfig}(undef, n)
    for i in 1:n
        if is_mpi && use_julia_amg
            # Avoid nesting HYPRE calls?
            amg_type = :aggregation
        else
            amg_type = :hypre
        end
        submodel = subsims[i].model
        linear_solver = select_reservoir_linear_solver(submodel, subdomain_precond)
        subconfigs[i] = simulator_config(subsims[i],
            max_timestep_cuts = inner_max_timestep_cuts,
            min_nonlinear_iterations = inner_min_nonlinear_iterations,
            tol_factor_final_iteration = inner_tol_final,
            check_before_solve = check_before_solve,
            linear_solver = linear_solver,
            info_level = Int(subdomain_info_level(i))
        )
        # Small hack for reservoir stuff
        set_default_cnv_mb!(subconfigs[i], submodel,
            tol_cnv = tol_multiplier(inner_tol_mul, :tol_cnv)*tol_cnv,
            tol_mb = tol_multiplier(inner_tol_mul, :tol_mb)*tol_mb,
            inc_tol_dp_abs = tol_multiplier(inner_tol_mul, :inc_tol_dp_abs)*inc_tol_dp_abs,
            inc_tol_dp_rel = tol_multiplier(inner_tol_mul, :inc_tol_dp_rel)*inc_tol_dp_rel,
            inc_tol_dz = tol_multiplier(inner_tol_mul, :inc_tol_dz)*inc_tol_dz,
            tol_dp_well = tol_multiplier(inner_tol_mul, :tol_dp_well)*tol_dp_well,
            tol_cnv_well = tol_multiplier(inner_tol_mul, :tol_cnv_well)*tol_cnv_well,
            tol_mb_well = tol_multiplier(inner_tol_mul, :tol_mb_well)*tol_mb_well
        )
        subconfigs[i][:id] = "Ω_$i"
        if !ismissing(inner_max_nonlinear_iterations) && isnothing(sim.storage[:subdomain_nldd_simulators][i])
            subconfigs[i][:max_nonlinear_iterations] = inner_max_nonlinear_iterations
        end
        add_option!(subconfigs[i], :subdomain_accept_unconverged, inner_accept_unconverged,
            "Keep the solution of a local solve that did not converge within the maximum number of iterations (no time-step cuts). The global stage then uses Newton rows linearized at this solution for the subdomain.", types = Bool)
        inner = sim.storage[:subdomain_nldd_simulators][i]
        if !isnothing(inner)
            # NLDD as solver for this subdomain. Global convergence must be
            # checked on the assembled system so that the subdomain Jacobian is
            # up to date for ASPEN.
            inner_cfg = simulator_config(inner;
                method = inner_method,
                gauss_seidel = true,
                info_level = -1,
                min_nonlinear_iterations = inner_min_nonlinear_iterations,
                linear_solver = linear_solver,
                tol_cnv = tol_multiplier(inner_tol_mul, :tol_cnv)*tol_cnv,
                tol_mb = tol_multiplier(inner_tol_mul, :tol_mb)*tol_mb,
                tol_cnv_well = tol_multiplier(inner_tol_mul, :tol_cnv_well)*tol_cnv_well,
                tol_mb_well = tol_multiplier(inner_tol_mul, :tol_mb_well)*tol_mb_well,
                tol_dp_well = tol_multiplier(inner_tol_mul, :tol_dp_well)*tol_dp_well,
                inc_tol_dp_abs = tol_multiplier(inner_tol_mul, :inc_tol_dp_abs)*inc_tol_dp_abs,
                inc_tol_dp_rel = tol_multiplier(inner_tol_mul, :inc_tol_dp_rel)*inc_tol_dp_rel,
                inc_tol_dz = tol_multiplier(inner_tol_mul, :inc_tol_dz)*inc_tol_dz,
                subdomain_tol_sufficient = false
            )
            if ismissing(inner_max_nonlinear_iterations)
                inner_cfg[:max_nonlinear_iterations] = subconfigs[i][:max_nonlinear_iterations]
            else
                inner_cfg[:max_nonlinear_iterations] = inner_max_nonlinear_iterations
            end
            add_option!(subconfigs[i], :subdomain_nldd, (inner, inner_cfg), "NLDD simulator and config used to solve this subdomain")
        end
    end
    add_option!(cfg, :config_subdomains, subconfigs, "Configs for subdomains")
    Jutul.overwrite_by_kwargs(cfg; check_before_solve = check_before_solve, kwarg...)
    for i in 1:n
        for k in [:failure_cuts_timestep, :relaxation]
            subconfigs[i][k] = cfg[k]
        end
        if haskey(subconfigs[i], :subdomain_nldd)
            # NLDD as subdomain solver: Same settings for its global and local solves
            inner_cfg = last(subconfigs[i][:subdomain_nldd])
            for c in vcat(inner_cfg, inner_cfg[:config_subdomains])
                for k in [:failure_cuts_timestep, :relaxation]
                    c[k] = cfg[k]
                end
            end
        end
    end
    cfg[:end_report] = cfg[:info_level] > -1
    return cfg
end

function partition_uniform_1d(nc, N)
    block_size = nc/N
    p = zeros(Integer, nc)
    prev = 0
    for i = 1:N
        if i == N
            next = nc
        else
            next = Integer(round(prev + block_size))
        end
        p[(prev+1):next] .= i
        prev = next
    end
    return p
end

function count_nldd_solve_statistics(reports)
    """
    Count NLDD solve statistics from simulation reports.
    Returns a named tuple with statistics or missing if required data is not available.
    """
    # Check if we have the basic structure we need
    if isempty(reports)
        return missing
    end
    
    # Initialize counters
    failure_count = 0
    nldd_count = 0
    total_count = 0
    subdomain_skipped = 0
    subdomain_solves = 0
    subdomain_total = 0
    
    # Track if we found any relevant data
    has_solve_data = false
    has_failure_data = false
    has_status_data = false
    
    for r in reports
        for m in r[:ministeps]
            if haskey(m, :steps)
                for mr in m[:steps]
                    # Count NLDD method usage
                    if haskey(mr, :converged) && haskey(mr, :local_solves_active)
                        has_solve_data = true
                        if !mr[:converged]
                            total_count += 1
                            nldd_count += mr[:local_solves_active]
                        end
                    end
                    
                    # Count failures and status
                    if haskey(mr, :subdomains)
                        num_sweeps = length(mr[:subdomains])
                        for sno in 1:num_sweeps
                            # Count subdomain failures
                            if haskey(mr, :subdomain_failures)
                                failures = mr[:subdomain_failures][sno]
                                if !isnothing(failures)
                                    has_failure_data = true
                                    failure_count += length(failures)
                                end
                            end
                            
                            # Count solve status
                            if haskey(mr, :solve_status)
                                has_status_data = true
                                for status in mr[:solve_status][sno]
                                    subdomain_skipped += status == local_solve_skipped
                                    subdomain_solves += status != local_already_converged && status != local_solve_skipped
                                    subdomain_total += 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    # Return missing if we don't have any meaningful data
    if !has_solve_data && !has_failure_data && !has_status_data
        return missing
    end
    
    subdomain_already_conv = subdomain_total - subdomain_solves - subdomain_skipped
    
    return (
        nldd_count = nldd_count,
        total_count = total_count,
        failure_count = failure_count,
        subdomain_skipped = subdomain_skipped,
        subdomain_solves = subdomain_solves,
        subdomain_total = subdomain_total,
        subdomain_already_conv = subdomain_already_conv,
        has_solve_data = has_solve_data,
        has_failure_data = has_failure_data,
        has_status_data = has_status_data
    )
end

function print_nldd_solve_statistics(stats, n_subdomains, method, info_level)
    """
    Print NLDD solve statistics if available and info_level allows it.
    """
    if ismissing(stats)
        return
    end
    
    if info_level > -1
        # Print NLDD method usage statistics
        if stats.has_solve_data
            Jutul.jutul_message("NLDD", "$(stats.nldd_count)/$(stats.total_count) solves used $method.")
        end
        
        # Print subdomain status statistics
        if stats.has_status_data
            Jutul.jutul_message("NLDD", "Subdomain status:\
            \n\t$(stats.subdomain_skipped)/$(stats.subdomain_total) local problems were skipped.\
            \n\t$(stats.subdomain_already_conv)/$(stats.subdomain_total) local problems were already converged.\
            \n\t$(stats.subdomain_solves)/$(stats.subdomain_total) local problems were solved.")
        end
    end
    
    # Always print failure information if there are failures
    if stats.has_failure_data && stats.failure_count > 0
        @info "$(stats.failure_count) subdomain solves failed out of $(stats.subdomain_solves) total local solves."
    end
end

function final_simulation_message(sim::NLDDSimulator, p, rec, t_elapsed, reports, timesteps, config, start_date, aborted)
    info_level = config[:info_level]
    if info_level > -1
        stats = subdomain_report_stats(sim, reports)
        if !ismissing(stats)
            Jutul.jutul_message("Subdomain", "statistics:")
            Jutul.print_iterations(stats, title = "Subdomains (Total)")
            Jutul.print_timing(stats, title = "Subdomain stats (Total)")
        end
    end
    
    # Count and print NLDD solve statistics
    nldd_stats = count_nldd_solve_statistics(reports)
    n_subdomains = length(sim.subdomain_simulators)
    print_nldd_solve_statistics(nldd_stats, n_subdomains, config[:method], info_level)
    
    final_simulation_message(sim.simulator, p, rec, t_elapsed, reports, timesteps, config, start_date, aborted)
end

function subdomain_report_stats(sim, reports)
    N = length(sim.subdomain_simulators)
    function get_cells(sim)
        if isa(sim.model, Jutul.MultiModel)
            d = sim.model.models.Reservoir.domain
        else
            d = sim.model.domain
        end
        return number_of_cells(d)
    end
    stats = Jutul.initialize_report_stats(reports)
    for outer_rep in reports
        # total_time += outer_rep[:total_time]
        for mini_rep in outer_rep[:ministeps]
            if haskey(mini_rep, :finalize_time)
                # total_finalize += mini_rep[:finalize_time]
            end
            if !haskey(mini_rep, :steps)
                return missing
            end
            for rep in mini_rep[:steps]
                if haskey(rep, :subdomains)
                    stats[:time] += sum(rep[:time_subdomains])
                    num_sweeps = length(rep[:subdomains])
                    for sno in 1:num_sweeps
                        sweep_rep = rep[:subdomains][sno]
                        if isnothing(sweep_rep)
                            continue
                        end
                        for i = 1:N
                            R = sweep_rep[i]
                            if isnothing(R)
                                continue
                            end
                            for step_rep in R
                                Jutul.ministep_report_stats!(stats, step_rep)
                            end
                        end
                    end
                end
            end
        end
    end
    Jutul.update_other_time_report_stats!(stats)
    return Jutul.output_report_stats(stats)
end

import JutulDarcy: reservoir_partition

function build_submodels(model, partition; active_global = missing, specialize = false)
    np = number_of_subdomains(partition)
    submodels = Vector{Any}(undef, np)
    p = Progress(np, desc = "Building $np submodels ")
    rmodel = reservoir_model(model)
    grid = physical_representation(rmodel.domain)
    N = get_neighborship(grid)
    fpos = get_facepos(N)

    for i = 1:np
        sm = submodel(model, partition, i, buffer = 1, minbatch = 1_000_000_000, active_global = active_global, facepos = fpos)
        if specialize
            sm = convert_to_immutable_storage(sm)
        end
        submodels[i] = sm
        next!(p)
    end
    return submodels
end

export bench_dd

function bench_dd(name, method = :fi; 
                                block_backend = true, 
                                steps = :full, info_level = 0, 
                                extra_timing = false, 
                                initial_dt = si_unit(:day), 
                                target_its = 8,
                                target_ds = Inf,
                                rtol = nothing,
                                wells = :ms,
                                local_rtol = rtol,
                                min_local_iterations = 0,
                                case = nothing,
                                mrst_data = missing,
                                global_linsolve = :cpr,
                                global_amg = JutulDarcy.default_amg_symbol(),
                                local_linsolve = :cpr,
                                global_krylov = :gmres,
                                local_krylov = :gmres,
                                inner_tol_mul = 1.0,
                                inner_tol_final = 1.0,
                                precompile = true,
                                tol_cnv = 1e-3,
                                tol_mb = 1e-7,
                                tol_dp_well = 1e-3,
                                tol_cnv_well = 10*tol_cnv,
                                tol_mb_well = 1e4*tol_mb,
                                inc_tol_dp_abs = Inf,
                                inc_tol_dp_rel = Inf,
                                inc_tol_dz = Inf,
                                timestepping = :iteration,
                                relaxation = NoRelaxation(),
                                maxstep = Inf,
                                restart = nothing,
                                ds_max = 0.2,
                                dp_max_abs = nothing,
                                dp_max_rel = 0.2,
                                p = nothing,
                                N = nothing,
                                kwarg...)
    do_print = info_level >= 0
    if do_print
        @info "Reading $name..."
    end
    if isnothing(case)
        if isfile(name)
            pth = name
        else
            pth = JutulDarcy.get_mrst_input_path(name)
        end
        base_path, filename = splitdir(pth)
        name, = splitext(filename)
        case, mrst_data = setup_case_from_mrst(pth,
                                        facility_grouping = :perwell,
                                        ds_max = ds_max,
                                        dp_max_rel = dp_max_rel,
                                        dp_max_abs = dp_max_abs,
                                        wells = wells,
                                        block_backend = block_backend)
    else
        base_path = ""
        case = deepcopy(case)
    end
    if do_print
        @info "Setting up $name..."
    end
    (; model, state0, parameters, forces) = case
    timesteps = case.dt
    if steps == :full
        dt = timesteps
    elseif steps == :first
        dt = timesteps[[1]];
    elseif steps == :mini
        dt = 0.001*timesteps[[1]]
    elseif steps isa Int64
        dt = timesteps[[steps]]
    elseif steps isa Vector{Int64} || steps isa UnitRange
        dt = timesteps[steps]
        if length(forces) == length(timesteps) && length(forces) > 1
            forces = forces[steps]
        end
    elseif steps isa Vector{Float64}
        dt = steps
    elseif steps isa Float64
        dt = [steps]
    else
        error()
    end
    ## Set up linear solver and preconditioner
    # Simulate
    lsolve = select_reservoir_linear_solver(model, global_linsolve, rtol = rtol, solver = global_krylov, amg_type = global_amg)
    function local_config(sim; kwarg...)
        tol_arg = (
            tol_cnv = tol_cnv,
            tol_mb = tol_mb,
            tol_cnv_well = tol_cnv_well,
            tol_mb_well = tol_mb_well,
            inc_tol_dp_abs = inc_tol_dp_abs,
            inc_tol_dp_rel = inc_tol_dp_rel,
            inc_tol_dz = inc_tol_dz,
            tol_dp_well = tol_dp_well
        )
        if method == :fi
            extra = NamedTuple()
        else
            extra = (inner_tol_mul = inner_tol_mul, inner_tol_final = inner_tol_final, tol_arg...)
        end
        cfg = simulator_config(sim, info_level = info_level,
                                extra_timing = extra_timing,
                                relaxation = relaxation,
                                linear_solver = lsolve; extra..., kwarg...)
        if method == :fi
            set_default_cnv_mb!(cfg, sim.model;
                tol_arg...
            )
        end
        return cfg
    end

    if method == :fi
        sim = Simulator(model, state0 = state0, parameters = deepcopy(parameters))
        cfg = local_config(sim; kwarg...)
    elseif method == :nldd || method == :aspen || method == :local_only
        nc = number_of_cells(model[:Reservoir].domain)
        if isnothing(p)
            p_name = "partition_$name.txt"
            # If N was specified we call the partioner first
            p = partition_from_N(model, parameters, N)
            # Check if it was provided
            if isnothing(p) && !ismissing(mrst_data)
                p = partition_from_mrst(model, mrst_data, do_print)
            end
            # Same folder as .mat file
            if isnothing(p)
                p_path = joinpath(base_path, p_name)
                p = partition_from_file(model, p_path, do_print)
            end
            # Inside NLDD path folder
            # if isnothing(p)
            #     p_path = joinpath(nldd_data_path, p_name)
            #     p = partition_from_file(model, p_path, do_print)
            # end
            # Fall back
            if isnothing(p)
                N = Int64(max(ceil(nc/1000), 2))
                if do_print
                    @info "Calling Metis with $N blocks since I did not find a partition anywhere and $N was not specified"
                end
                p = partition_internal_metis(model, parameters, N, do_print)
            end
        else
            if do_print
                @info "Got partition p as kwarg"
            end
            @assert length(p) == nc
            p::AbstractVecOrMat
            p = vec(p)
        end
        N = maximum(p)
        if do_print
            @info "Partition ready with $N coarse blocks."
        end
        @assert length(p) == nc "Expected partition with $nc entries, was $(length(p))"
        @assert maximum(p) == length(unique(p))
        if N == 1
            @warn "Partition has one entry for all cells. Possible mistake?"
        end
        mpart = reservoir_partition(model, p);
        submodels = build_submodels(model, mpart);
        sim = NLDDSimulator(model, mpart, submodels = submodels, state0 = state0, parameters = deepcopy(parameters));
        cfg = local_config(sim; method = method, kwarg...)
        for c in cfg[:config_subdomains]
            c[:linear_solver] = select_reservoir_linear_solver(model, local_linsolve, rtol = local_rtol, v = 0, solver = local_krylov)
            c[:min_nonlinear_iterations] = min_local_iterations
            # c[:max_nonlinear_iterations] = 20
            c[:relaxation] = cfg[:relaxation]
        end
    else
        error("Method $method is not supported.")
    end
    if timestepping == :iteration
        t_base = TimestepSelector(initial_absolute = initial_dt, max = maxstep)
        t_its = IterationTimestepSelector(target_its)
        if isfinite(target_ds)
            t_sat = VariableChangeTimestepSelector(:Saturations, target_ds, relative = false, reduction = :max, model = :Reservoir)
            tstepper = [t_base, t_its, t_sat]
        else
            tstepper = [t_base, t_its]
        end
    else
        @assert isnothing(timestepping)
        t_base = TimestepSelector(max = maxstep)
        tstepper = [t_base]
    end
    cfg[:timestep_selectors] = tstepper
    opth = cfg[:output_path]
    if do_print
        @info "Simulating $name with $method..."
    end
    if precompile
        minicase = case[1:1]
        minicase.dt[1] *= 0.01
        cfg[:output_path] = nothing
        old_il = cfg[:info_level]
        old_er = cfg[:end_report]
        # Disable printing for precompilation
        cfg[:info_level] = -1
        cfg[:end_report] = false

        simulate!(sim, minicase.dt,
            forces = minicase.forces,
            config = cfg,
            restart = nothing,
            state0 = state0
        )
        cfg[:info_level] = old_il
        cfg[:end_report] = old_er
        cfg[:output_path] = opth
    end
    result = simulate!(sim, dt,
        forces = forces,
        config = cfg,
        restart = restart,
        state0 = state0)
    return ReservoirSimResult(model, result, forces, config = cfg, case = case, sim = sim, mrst = mrst_data, partition = p)
end

export nldd_output_path
function nldd_output_path(casename, extra = ""; method = "nldd")
    return jutul_output_path(method, subfolder = joinpath("nldd", extra, casename))
end

function partition_from_N(model, parameters, N)
    if isnothing(N)
        p = nothing
    else
        p = partition_internal_metis(model, parameters, N)
        jutul_message("NLDD-partitioner", "N = $N provided, calling Metis.")
    end
    return p
end

function partition_from_mrst(model, mrst_data, do_print = true)
    p = nothing
    if haskey(mrst_data, "extra") && length(mrst_data["extra"]) > 0
        e = first(mrst_data["extra"])
        if do_print
            jutul_message("NLDD-partitioner", "Partition found in MRST export.")
        end
        if haskey(e, "partition")
            p = vec(Int64.(e["partition"]))
            if length(p) == 0
                p = nothing
            else
                if do_print
                    jutul_message("NLDD-partitioner", "Partition found in MRST export (under extra).")
                end
            end
        end
    elseif haskey(mrst_data, "partition")
        p = vec(Int64.(mrst_data["partition"]))
        if do_print
            jutul_message("NLDD-partitioner", "Partition found in MRST export.")
        end
    else
        p = nothing
    end
    return p
end

function partition_from_file(model, p_path, do_print = true)
    if isfile(p_path)
        p_raw = readlines(p_path)
        p = parse.(Int64, p_raw)
        if minimum(p) == 0
            @. p += 1
        end
        if do_print
            @info "Partition found in file $p_path"
        end
    else
        p = nothing
    end
    return p
end

function bench_dd(name, methods::Vector{Symbol}, kwarg...)
    for m in methods
        bench_dd(name, m; kwarg...)
    end
end

export aspen_forces
function aspen_forces(forces::Vector, submodels)
    return map(f -> aspen_forces(f, submodels), forces)
end

function aspen_forces(forces, submodels)
    return (outer = forces, subdomains = map(x -> subforces(forces, x), submodels))
end

function Jutul.progress_recorder(sim::NLDDSimulator)
    return sim.simulator.storage.recorder
end

function set_default_cnv_mb!(config::JutulConfig, sim::NLDDSimulator; kwarg...)
    set_default_cnv_mb!(config, sim.simulator.model; kwarg...)
    sub_kwarg = scale_subdomain_tolerances(config[:inner_tol_mul]; kwarg...)
    for i in eachindex(sim.subdomain_simulators)
        cfg = config[:config_subdomains][i]
        m = sim.subdomain_simulators[i].model
        set_default_cnv_mb!(cfg, m; sub_kwarg...)
        if haskey(cfg, :subdomain_nldd)
            inner_sim, inner_cfg = cfg[:subdomain_nldd]
            set_default_cnv_mb!(inner_cfg, inner_sim; sub_kwarg...)
        end
    end
end

"""
    local_well_reservoir_cells(model)

Reservoir cells perforated by each well, in the cell numbering of `model` (which
may be the model of a subdomain or an MPI process), taken from the well cross
terms.
"""
function local_well_reservoir_cells(model::MultiModel)
    out = OrderedDict{Symbol, Vector{Int}}()
    for ctp in model.cross_terms
        ct = ctp.cross_term
        if ct isa JutulDarcy.ReservoirFromWellFlowCT
            if ctp.target == :Reservoir
                w = ctp.source
            else
                w = ctp.target
            end
            out[w] = unique(vcat(get(out, w, Int[]), collect(ct.reservoir_cells)))
        end
    end
    return out
end

local_well_reservoir_cells(model) = OrderedDict{Symbol, Vector{Int}}()

"""
    local_reservoir_partition(model, nblocks; active = missing)

Partition the reservoir of `model` into `nblocks` blocks with Metis, using the
connections in the cell numbering of the model itself. This works for the
models of subdomains and MPI processes, where `partition_reservoir` (which uses
the connections of the global mesh) does not. Cells that are not unknowns
(`active[c] = false`, e.g. ghost or boundary cells) are left out of the graph
and assigned to the block of a neighbor. All cells perforated by a well end up
in the same block.
"""
function local_reservoir_partition(model, nblocks::Int; active = missing, parameters = missing, conn = :unit)
    rmodel = JutulDarcy.reservoir_model(model)
    D = rmodel.domain
    nc = number_of_cells(D)
    if ismissing(active)
        active = trues(nc)
    end
    wells = local_well_reservoir_cells(model)
    if nblocks <= 1
        p = ones(Int, nc)
    else
        N = get_neighborship(physical_representation(D))
        ia = findall(active)
        local_index = zeros(Int, nc)
        local_index[ia] = eachindex(ia)
        keep = [active[N[1, f]] && active[N[2, f]] for f in axes(N, 2)]
        N_a = local_index[N[:, keep]]
        edge_weights = local_partition_weights(parameters, conn, size(N, 2))[keep]
        groups = Vector{Vector{Int}}()
        for cells in values(wells)
            g = unique(local_index[filter(c -> active[c], cells)])
            length(g) > 1 && push!(groups, g)
        end
        p_a = Jutul.partition_hypergraph(N_a, nblocks, Jutul.MetisPartitioner(); num_nodes = length(ia), groups = groups, edge_weights = edge_weights)
        p = zeros(Int, nc)
        p[ia] = p_a
        # Inactive cells get the block of a neighbor
        for f in axes(N, 2)
            l, r = N[1, f], N[2, f]
            if p[l] == 0 && p[r] > 0
                p[l] = p[r]
            elseif p[r] == 0 && p[l] > 0
                p[r] = p[l]
            end
        end
        p[p .== 0] .= 1
        # Keep the cells of each well in a single block
        for cells in values(wells)
            blocks = unique(p[cells])
            for b in blocks[2:end]
                p[p .== b] .= blocks[1]
            end
        end
        p = Jutul.compress_partition(p)
    end
    return local_multimodel_partition(model, p, wells)
end

# Integer connection weights for Metis, as in partitioner_input (:unit or :trans)
function local_partition_weights(parameters, conn, nf)
    if conn == :unit || ismissing(parameters)
        return ones(Int, nf)
    end
    conn == :trans || error("Local partition supports conn = :unit or :trans, was $conn")
    prm = haskey(parameters, :Reservoir) ? parameters[:Reservoir] : parameters
    T = max.(Jutul.value.(prm[:Transmissibilities]), 1e-20)
    length(T) == nf || error("Expected $nf transmissibilities, found $(length(T)).")
    T = length(T)*T./sum(T)
    return Int.(ceil.(T))
end

function local_multimodel_partition(model::MultiModel, p, wells)
    part = Dict{Symbol, Any}()
    for key in keys(model.models)
        if key == :Reservoir
            part[key] = Jutul.SimplePartition(p)
        elseif haskey(wells, key)
            blocks = unique(p[wells[key]])
            @assert length(blocks) == 1 "All cells of well $key must be in the same block, found $blocks."
            part[key] = only(blocks)
        else
            # Well controller of a split well
            w = Symbol(string(key)[1:end-5])
            part[key] = part_for_well(p, wells, w, key)
        end
    end
    return Jutul.SimpleMultiModelPartition(part, :Reservoir)
end

function part_for_well(p, wells, w, key)
    haskey(wells, w) || error("Could not find the well of model $key for partitioning.")
    return only(unique(p[wells[w]]))
end

local_multimodel_partition(model, p, wells) = Jutul.SimplePartition(p)

"""
    subdomain_nldd_simulator(subsim, sub_model, state0, parameters, nblocks, parent_model; conn = :unit)

NLDD simulator with `nblocks` blocks that solves the subdomain of `subsim`. The
boundary cells of the subdomain are kept fixed, and the NLDD simulator shares
`subsim` (state, linear system) as its outer simulator. `parent_model` is the
model the subdomain was taken from. `conn` gives the edge weights of the
partition (see `local_reservoir_partition`).
"""
function subdomain_nldd_simulator(subsim, sub_model, state0, parameters, nblocks, parent_model; conn = :unit)
    gmap = JutulDarcy.reservoir_model(sub_model).domain.global_map
    active = .!gmap.cell_is_boundary
    partition = local_reservoir_partition(sub_model, nblocks; active = active, parameters = parameters, conn = conn)
    blocks = inner_block_submodels(parent_model, sub_model, partition)
    case = JutulCase(sub_model, parameters = parameters, state0 = state0)
    return NLDDSimulator(case, partition;
        submodels = blocks,
        outer_simulator = subsim,
        setup_aspen = false
    )
end

"""
    inner_block_submodels(parent_model, sub_model, partition)

Models of the blocks in `partition` (given in the cell numbering of `sub_model`).
Submodels of a subdomain model cannot be built directly, so the blocks are
taken from `parent_model` (which `sub_model` is a subdomain of), and the
reservoir global map is then changed to point to the cells of `sub_model`.
"""
function inner_block_submodels(parent_model, sub_model, partition)
    gmap = JutulDarcy.reservoir_model(sub_model).domain.global_map
    active = .!gmap.cell_is_boundary
    p_sub = Jutul.main_partition(partition).partition
    K = maximum(p_sub)
    rparent = JutulDarcy.reservoir_model(parent_model)
    nc_parent = number_of_cells(rparent.domain)
    # Cells outside the subdomain (and models not in it) go to an extra block
    # that is never built.
    p_parent = fill(K + 1, nc_parent)
    for c in eachindex(p_sub)
        if active[c]
            p_parent[gmap.cells[c]] = p_sub[c]
        end
    end
    active_parent = p_parent .<= K
    if parent_model isa MultiModel
        part = Dict{Symbol, Any}()
        for key in keys(parent_model.models)
            if key == :Reservoir
                part[key] = Jutul.SimplePartition(p_parent)
            elseif haskey(partition.partition, key)
                part[key] = partition.partition[key]
            else
                part[key] = K + 1
            end
        end
        parent_partition = Jutul.SimpleMultiModelPartition(part, :Reservoir)
    else
        parent_partition = Jutul.SimplePartition(p_parent)
    end
    N = get_neighborship(physical_representation(rparent.domain))
    fpos = get_facepos(N)
    cell_index = Dict(c => i for (i, c) in enumerate(gmap.cells))
    face_index = Dict(f => i for (i, f) in enumerate(gmap.faces))
    blocks = Vector{Any}(undef, K)
    for j in 1:K
        sm = submodel(parent_model, parent_partition, j, buffer = 1, minbatch = 1_000_000_000,
            active_global = active_parent, facepos = fpos)
        blocks[j] = reservoir_global_map_to_local(sm, cell_index, face_index)
    end
    return blocks
end

function reservoir_global_map_to_local(model::MultiModel, cell_index, face_index)
    model.models[:Reservoir] = reservoir_global_map_to_local(model.models[:Reservoir], cell_index, face_index)
    return model
end

function reservoir_global_map_to_local(model::SimulationModel, cell_index, face_index)
    D = model.domain
    M = D.global_map
    cells = [cell_index[c] for c in M.cells]
    faces = [get(face_index, f, 0) for f in M.faces]
    M_new = Jutul.FiniteVolumeGlobalMap(cells, faces, M.cell_is_boundary, variables_always_active = M.variables_always_active)
    D_new = Jutul.DiscretizedDomain(D.representation, D.discretizations, D.entities, M_new)
    T = typeof(model)
    params = collect(T.parameters)
    params[1] = typeof(D_new)
    fields = Any[getfield(model, f) for f in fieldnames(T)]
    fields[1] = D_new
    return SimulationModel{params...}(fields...)
end

"""
    tol_multiplier(mul, name)

Factor for the subdomain tolerance `name` (e.g. `:tol_cnv`). `mul` is either a
number used for all tolerances, or a NamedTuple with factors for some
tolerances (1.0 for the others).
"""
tol_multiplier(mul::Real, name::Symbol) = mul
tol_multiplier(mul::NamedTuple, name::Symbol) = Float64(get(mul, name, 1.0))

"""
    scale_subdomain_tolerances(mul; kwarg...)

Multiply all tolerance arguments (`tol_*` and `inc_tol_*`) by the factor given
by `inner_tol_mul` (see `tol_multiplier`).
"""
function scale_subdomain_tolerances(mul; kwarg...)
    out = Dict{Symbol, Any}()
    for (k, v) in kwarg
        s = string(k)
        if v isa Real && (startswith(s, "tol_") || startswith(s, "inc_tol_"))
            v = tol_multiplier(mul, k)*v
        end
        out[k] = v
    end
    return out
end

function Jutul.preprocess_forces(sim::NLDDSimulator, forces)
    per_step = forces isa Vector
    submodels = map(Jutul.get_simulator_model, sim.subdomain_simulators)
    forces = aspen_forces(forces, submodels)
    return (forces = forces, forces_per_step = per_step)
end

function Jutul.get_simulator_model(sim::NLDDSimulator)
    Jutul.get_simulator_model(sim.simulator)
end

function Jutul.get_simulator_storage(sim::NLDDSimulator)
    Jutul.get_simulator_storage(sim.simulator)
end
