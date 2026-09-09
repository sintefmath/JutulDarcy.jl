using Jutul, JutulDarcy, Test, LinearAlgebra
import Jutul.DictOptimization: finite_difference_gradient_entry, optimizer_devectorize

# Gradient checks for the well-control optimization layer, same structure as
# examples/workflow/control_optimization.jl (rate targets + BHP-cap limits on the
# injectors, BHP targets on the producers, NPV objective) on a hard-coarsened Egg
# (3x3x1, first 1.5 years, two control periods) for test speed. We check the
# period-2 gradient for six DOFs; with these settings the INJECT3 BHP cap partly
# binds in period 2 and the INJECT7 cap stays slack, so the checks cover an
# active and an inactive operating limit.
#
# di_sparse = true is the default. Its dR/dX sparsity pattern is detected once
# and frozen; wc_rebuild_case adds a value-neutral coupling term so a limit that
# is slack at detection still gets a correct gradient once it binds. The second
# testset guards that: the pattern is frozen at the base schedule (INJECT7
# slack), then re-evaluated with the injectors pushed to max rate (INJECT7 now
# capped).

function egg_wco_case(; dims = (3, 3, 1), years = 1.5, nperiods = 2)
    day = si_unit(:day)
    bar = si_unit(:bar)
    egg_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("EGG")
    fine = setup_case_from_data_file(joinpath(egg_dir, "EGG.DATA"))
    n5 = count(<=(years * si_unit(:year)), cumsum(fine.dt))
    egg = coarsen_reservoir_case(fine[1:n5], dims, method = :ijk)
    nstep = length(egg.dt)

    t = cumsum(egg.dt)
    step_to_control = ceil.(Int, nperiods .* t ./ t[end])
    edges = [0; findall(diff(step_to_control) .!= 0); nstep]
    control_periods = [(edges[i] + 1):edges[i + 1] for i in 1:nperiods]

    injectors = [Symbol("INJECT$i") for i in 1:8]
    producers = [Symbol("PROD$i") for i in 1:4]
    controls = NamedTuple[]
    for w in injectors
        push!(controls, (well = w, quantity = :rate, abs_min = 0.0 / day, abs_max = 150.0 / day))
        push!(controls, (well = w, quantity = :bhp,  abs_min = 400 * bar,  abs_max = 420 * bar))
    end
    for w in producers
        push!(controls, (well = w, quantity = :bhp, abs_min = 380 * bar, abs_max = 400 * bar))
    end
    npv_arg = (oil_price = 50.0, water_price = -2.0, water_cost = 2.0, oil_cost = 0.0,
               discount_rate = 0.1, scale = 1e8)
    copt = setup_well_control_optimization(egg, control_periods, controls;
        objective = :npv, npv_arg = npv_arg, verbose = false)
    return (; copt, egg, control_periods, injectors, producers, day, bar)
end

"Scalar index of DOF `(well, quantity)` for control period `period` in the flat control vector."
function dof_index(copt, well, quantity, period)
    off = 0
    for d in copt.dofs
        if d.well == well && d.quantity == quantity
            k = d.constant ? 1 : findfirst(==(period), d.periods)
            k === nothing && error("period $period not active for $well/$quantity")
            return off + k
        end
        off += d.constant ? 1 : length(d.periods)
    end
    error("no DOF for $well/$quantity")
end

"Central FD of the objective wrt scaled DOF `i`, from two `finite_difference_gradient_entry` calls."
fd_central(prob, x, i; h = 1e-5) =
    (finite_difference_gradient_entry(prob, x; index = i, eps = h) +
     finite_difference_gradient_entry(prob, x; index = i, eps = -h)) / 2

@testset "well control optimization" begin
    s = egg_wco_case()
    copt = s.copt
    period2_steps = s.control_periods[2]
    on_bhp(ws, w) = any(==(:bhp), ws.wells[w][:control][period2_steps])

    sim_arg = (rtol = 1e-5, tol_cnv = 1e-5)   # match the example (info_level = -1 is forced internally)

    @testset "DOF setup" begin
        @test length(copt.dofs) == 2 * length(s.injectors) + length(s.producers)
        for d in copt.dofs
            primary = d.well in s.injectors ? :rate : :bhp
            @test d.is_target == (d.quantity == primary)
            @test d.abs_min < d.abs_max
        end
        prm = copt.dopt.parameters
        f1 = copt(prm).forces[first(s.control_periods[1])][:Facility]
        @test f1.control[:INJECT3].target isa TotalRateTarget
        @test f1.control[:INJECT3].target.value ≈ prm["INJECT3"]["rate"][1]
        @test f1.limits[:INJECT3][:bhp] ≈ prm["INJECT3"]["bhp"][1]
        @test f1.control[:PROD1].target isa BottomHolePressureTarget
        @test f1.control[:PROD1].target.value ≈ prm["PROD1"]["bhp"][1]
    end

    prob   = JutulDarcy.well_control_optimization_problem(copt; simulator_arg = sim_arg)  # sparse (di_sparse = true)
    prob_d = JutulDarcy.well_control_optimization_problem(copt; simulator_arg = sim_arg,
        backend_arg = (di_sparse = false,))                                              # dense reference

    checks = [(:INJECT3, :rate), (:INJECT3, :bhp), (:INJECT7, :rate),
              (:INJECT7, :bhp), (:PROD1, :bhp), (:PROD3, :bhp)]
    x0 = copy(prob.x0)

    @testset "base schedule, period-2 gradients" begin
        f, g = prob(x0)
        @test isfinite(f) && f > 0
        _, gd = prob_d(x0)
        @test isapprox(g, gd; atol = 1e-8, rtol = 1e-5)          # sparse == dense
        gs = max(norm(gd, Inf), 1e-8)

        ws0, = simulate_reservoir(copt(copt.dopt.parameters); info_level = -1)
        @test on_bhp(ws0, :INJECT3)                              # BHP cap partly binds in period 2
        @test !on_bhp(ws0, :INJECT7)                             # BHP cap stays slack

        for (w, q) in checks
            i = dof_index(copt, w, q, 2)
            gfd = fd_central(prob_d, x0, i)
            inactive_limit = q == :bhp && w in s.injectors && !on_bhp(ws0, w)
            if inactive_limit
                @test abs(gd[i]) < 1e-4 * gs                     # inactive limit: ~0 adjoint
                @test abs(gfd)   < 1e-4 * gs                     #                 ~0 FD
            else
                near_zero = max(abs(gd[i]), abs(gfd)) < 1e-4 * gs
                @test near_zero || isapprox(gd[i], gfd; rtol = 1e-2, atol = 1e-4)
            end
        end
    end

    @testset "frozen sparse pattern stays correct when a limit activates" begin
        x_hi = copy(x0)
        for w in s.injectors, p in 1:length(s.control_periods)
            x_hi[dof_index(copt, w, :rate, p)] = 1.0             # push every injector to 150 m3/day
        end
        wsh, = simulate_reservoir(copt(optimizer_devectorize(prob, x_hi)); info_level = -1)
        @test on_bhp(wsh, :INJECT7)                              # its BHP cap now binds

        _, gh  = prob(x_hi)                                      # sparse: pattern was frozen at x0 (INJECT7 slack)
        _, gdh = prob_d(x_hi)
        @test isapprox(gh, gdh; atol = 1e-8, rtol = 1e-5)        # still equals dense

        i = dof_index(copt, :INJECT7, :bhp, 2)
        gfd = fd_central(prob_d, x_hi, i)
        @test isapprox(gh[i], gfd; rtol = 0.2, atol = 1e-4)      # ... and matches FD
        @test abs(gh[i]) > 1e-3 * max(norm(gdh, Inf), 1e-8)      # genuinely nonzero now
    end
end
