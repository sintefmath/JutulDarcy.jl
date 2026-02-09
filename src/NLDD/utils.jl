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
    add_option!(cfg, :method, method, "Method to use", values = [:nldd, :aspen], types = Symbol)
    add_option!(cfg, :debug_checks, false, "Enable extra expensive checks", types = Bool)
    add_option!(cfg, :solve_subdomains, true, "Solve subdomains", types = Bool)
    add_option!(cfg, :aspen_full_increment, false, "Solve full ASPEN update", types = Bool)
    add_option!(cfg, :strategy, DefaultNLDDStrategy(), "Strategy to use for applying NLDD/ASPEN")
    same_tol = inner_tol_final <= 1.0 && inner_tol_mul <= 1.0
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
            amg_type = :smoothed_aggregation
        else
            amg_type = JutulDarcy.default_amg_symbol()
        end
        submodel = subsims[i].model
        linear_solver = reservoir_linsolve(submodel, subdomain_precond)
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
            tol_cnv = inner_tol_mul*tol_cnv,
            tol_mb = inner_tol_mul*tol_mb,
            inc_tol_dp_abs = inner_tol_mul*inc_tol_dp_abs,
            inc_tol_dp_rel = inner_tol_mul*inc_tol_dp_rel,
            inc_tol_dz = inner_tol_mul*inc_tol_dz,
            tol_dp_well = inner_tol_mul*tol_dp_well,
            tol_cnv_well = inner_tol_mul*tol_cnv_well,
            tol_mb_well = inner_tol_mul*tol_mb_well
        )
        subconfigs[i][:id] = "Ω_$i"
    end
    add_option!(cfg, :config_subdomains, subconfigs, "Configs for subdomains")
    Jutul.overwrite_by_kwargs(cfg; check_before_solve = check_before_solve, kwarg...)
    for i in 1:n
        for k in [:failure_cuts_timestep, :relaxation]
            subconfigs[i][k] = cfg[k]
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
    n = length(sim.subdomain_simulators)
    failure_count = 0
    count = 0
    nldd_count = 0
    total_count = 0
    subdomain_skipped = 0
    subdomain_solves = 0
    subdomain_total = 0
    for r in reports
        for m in r[:ministeps]
            if haskey(m, :steps)
                for mr in m[:steps]
                    if haskey(mr, :subdomain_failures)
                        failures = mr[:subdomain_failures]
                        if isnothing(failures)
                            continue
                        end
                        failure_count += length(failures)
                        count += 1
                    end
                    if !mr[:converged]
                        total_count += 1
                        nldd_count += mr[:local_solves_active]
                    end
                    if haskey(mr, :solve_status)
                        for status in mr[:solve_status]
                            subdomain_skipped += status == local_solve_skipped
                            subdomain_solves += status != local_already_converged && status != local_solve_skipped
                            subdomain_total += 1
                        end
                    end
                end
            end
        end
    end
    if info_level > -1
        subdomain_already_conv = subdomain_total-subdomain_solves-subdomain_skipped
        Jutul.jutul_message("NLDD", "$nldd_count/$total_count solves used $(config[:method]).")
        Jutul.jutul_message("NLDD", "Subdomain status:\n\t$subdomain_skipped/$subdomain_total local solves were skipped.\n\t$subdomain_already_conv/$subdomain_total local solves were already converged.\n\t$subdomain_solves/$subdomain_total local solves were solved.")
    end
    if failure_count > 0
        @info "$failure_count subdomain solves failed out of $(n*count) total local solves."
    end
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
    lsolve = reservoir_linsolve(model, global_linsolve, rtol = rtol, solver = global_krylov, amg_type = global_amg)
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
            c[:linear_solver] = reservoir_linsolve(model, local_linsolve, rtol = local_rtol, v = 0, solver = local_krylov)
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
    for i in eachindex(sim.subdomain_simulators)
        cfg = config[:config_subdomains][i]
        m = sim.subdomain_simulators[i].model
        set_default_cnv_mb!(cfg, m; kwarg...)
    end
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
