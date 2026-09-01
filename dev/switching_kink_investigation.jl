# Investigation: the "switching kink" in d(NPV)/d(bhp limit)
# ==========================================================
#
# Background (see dev/well_control_optimization_demo.jl): for the demo case,
# the adjoint gradient w.r.t. the producer bhp *limit* in control period 1
# (dof 4) disagrees with central finite differences — the adjoint is negative
# while all FD slopes are positive, and one-sided FDs behave like there is a
# *discontinuity* (jump), not just a kink, in NPV(bhp_limit) at the base value.
#
# This script maps out NPV as a function of the period-1 bhp limit in detail:
#   1. Coarse sweep over the whole box to see the global shape.
#   2. Fine zoom around the base value (120 bar) to distinguish jump vs kink.
#   3. Adjoint gradients at selected points, compared against numerical slopes
#      of the sweep curve (piecewise, away from any jumps).
#   4. Control-regime diagnostics (which report steps operate on bhp vs orat)
#      to correlate every jump/kink with a switch-pattern change.
#
# Figures are written to dev/switching_kink_figs/*.png.
#
# Run:  julia --project=. dev/switching_kink_investigation.jl
# Needs branches: Jutul fix/adjoint-well-control-switching + this JutulDarcy branch.

using Jutul, JutulDarcy, Printf, LinearAlgebra, Statistics

const HAVE_MAKIE = !isnothing(Base.find_package("GLMakie"))
if HAVE_MAKIE
    using GLMakie
    GLMakie.activate!()
end

const bar = 1e5
const day = si_unit(:day)

# ---------------------------------------------------------------------------
# Identical case to well_control_optimization_demo.jl
# ---------------------------------------------------------------------------
function demo_case(; nstep = 24, dtdays = 20.0)
    mesh = CartesianMesh((15, 15, 2), (600.0, 600.0, 20.0))
    domain = reservoir_domain(mesh; permeability = 1e-13, porosity = 0.22)
    Inj  = setup_vertical_well(domain, 1, 1;   name = :INJ)
    Prod = setup_vertical_well(domain, 15, 15; name = :PROD)

    rhoS = [1000.0, 780.0]
    sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()); reference_densities = rhoS)
    model, parameters = setup_reservoir_model(domain, sys; wells = [Inj, Prod], extra_out = true)
    replace_variables!(model, PhaseMassDensities = ConstantCompressibilityDensities(
        p_ref = 200*bar, density_ref = rhoS, compressibility = [1e-5/bar, 8e-5/bar]))

    state0 = setup_reservoir_state(model; Pressure = 200*bar, Saturations = [0.1, 0.9])
    dt = repeat([dtdays*day], nstep)

    base_rate = 900.0/day
    controls = Dict(
        :INJ  => InjectorControl(TotalRateTarget(base_rate), [1.0, 0.0], density = rhoS[1]),
        :PROD => ProducerControl(SurfaceOilRateTarget(-600.0/day)),
    )
    limits = Dict(
        :INJ  => (bhp = 350*bar,),
        :PROD => (bhp = 120*bar,),
    )
    forces = setup_reservoir_forces(model; control = controls, limits = limits)
    return JutulCase(model, dt, forces; state0 = state0, parameters = parameters)
end

case = demo_case()
periods = [1:8, 9:16, 17:24]

# Single DOF: the PROD bhp limit in period 1 only. Periods 2-3 keep the base
# 120 bar limit. scaler = missing -> the optimization variable is the physical
# bhp value in Pa, so gradients need no rescaling.
b_lo, b_hi = 95*bar, 185*bar
copt = setup_well_control_optimization(case, periods,
    [(well = :PROD, quantity = :bhp, abs_min = b_lo, abs_max = b_hi,
      periods = [1], scaler = missing)];
    objective = :npv,
    npv_arg = (oil_price = 60.0, water_price = -6.0, water_cost = 3.0, discount_rate = 0.1),
    verbose = false)

prob = JutulDarcy.well_control_optimization_problem(copt)

npv(b)      = first(prob([b]; gradient = false))         # forward only (substates)
npv_grad(b) = prob([b])                                   # (f, [dNPV/db]) adjoint

# Fresh forward sim for control-pattern diagnostics (report-step resolution)
function control_diag(b)
    prm = deepcopy(copt.dopt.parameters)
    prm["PROD"]["bhp"] = [b]
    c = copt.setup_function(prm, missing)
    r = simulate_reservoir(c; info_level = -1)
    ctrl = r.wells[:PROD][:control]
    return (
        ctrl = ctrl,
        n_bhp_p1 = count(==(:bhp), ctrl[1:8]),
        n_bhp = count(==(:bhp), ctrl),
        bhp = r.wells[:PROD][:bhp],
    )
end

fmtbar(b) = @sprintf("%.3f", b/bar)

# ---------------------------------------------------------------------------
# 1) Coarse sweep over the box
# ---------------------------------------------------------------------------
println("== 1) Coarse sweep ==")
bs_coarse = collect(range(b_lo, b_hi; step = 2.5*bar))
f_coarse = Float64[]
for (i, b) in enumerate(bs_coarse)
    push!(f_coarse, npv(b))
    @printf("  %2d/%2d  bhp = %7s bar   NPV = %.7e\n", i, length(bs_coarse), fmtbar(b), f_coarse[end])
    flush(stdout)
end

# ---------------------------------------------------------------------------
# 2) Fine zoom around the base value
# ---------------------------------------------------------------------------
println("== 2) Zoom sweep around 120 bar ==")
bs_zoom = collect(range(118*bar, 122*bar; length = 81))
f_zoom = [npv(b) for b in bs_zoom]

# Detect jump candidates: adjacent differences far above the local median
df = diff(f_zoom)
db = step(range(118*bar, 122*bar; length = 81))
med = median(abs.(df))
jump_ix = findall(d -> abs(d) > 8*med, df)
println("  median |dNPV| per $(fmtbar(db)) bar: $(med)")
for j in jump_ix
    @printf("  JUMP candidate between %s and %s bar: dNPV = %+.4e (%.0fx median)\n",
        fmtbar(bs_zoom[j]), fmtbar(bs_zoom[j+1]), df[j], abs(df[j])/med)
end

# ---------------------------------------------------------------------------
# 3) Adjoint vs numerical slopes
# ---------------------------------------------------------------------------
println("== 3) Adjoint gradients ==")
bs_adj = vcat(collect(range(100*bar, 180*bar; step = 10*bar)), [119.0*bar, 120.0*bar, 121.0*bar])
sort!(unique!(bs_adj))
g_adj = Float64[]
f_adj = Float64[]
for b in bs_adj
    f, g = npv_grad(b)
    push!(f_adj, f); push!(g_adj, g[1])
    @printf("  bhp = %7s bar   NPV = %.6e   adjoint dNPV/db = %+.6e (per Pa)\n",
        fmtbar(b), f, g[1])
    flush(stdout)
end

# Numerical slopes from the zoom grid, one-sided at each adjoint point near 120
println("== One-sided slopes at 120 bar (from dedicated evals) ==")
b0 = 120.0*bar
f0 = npv(b0)
println("  f(120 bar) = $f0")
for h in (0.05, 0.1, 0.2, 0.5, 1.0) .* bar
    fp = npv(b0 + h); fm = npv(b0 - h)
    @printf("  h = %5.2f bar  fwd (f(x+h)-f)/h = %+.5e   bwd (f-f(x-h))/h = %+.5e   jump est f-f(x-h) = %+.4e\n",
        h/bar, (fp - f0)/h, (f0 - fm)/h, f0 - fm)
end
i120 = argmin(abs.(bs_adj .- b0))
@printf("  adjoint at 120 bar: %+.5e per Pa\n", g_adj[i120])

# ---------------------------------------------------------------------------
# 4) Control-pattern diagnostics
# ---------------------------------------------------------------------------
println("== 4) Control patterns ==")
diag_bs = sort(unique(vcat(collect(range(100*bar, 180*bar; step = 10*bar)),
    [119.0, 119.8, 119.9, 119.95, 120.0, 120.05, 120.1, 120.5, 121.0] .* bar,
    [bs_zoom[j] for j in jump_ix], [bs_zoom[j+1] for j in jump_ix])))
diag = Dict{Float64, Any}()
for b in diag_bs
    d = control_diag(b)
    diag[b] = d
    @printf("  bhp = %8s bar   steps on :bhp (period1/total) = %d/%d   ctrl = %s\n",
        fmtbar(b), d.n_bhp_p1, d.n_bhp, join(first.(string.(d.ctrl)), ""))
    flush(stdout)
end

# ---------------------------------------------------------------------------
# 5) Plots
# ---------------------------------------------------------------------------
if HAVE_MAKIE
    figdir = joinpath(@__DIR__, "switching_kink_figs")
    mkpath(figdir)
    MUSD = 1e6

    # Fig 1: global NPV curve + adjoint tangents + control-regime coloring
    fig1 = Figure(size = (1100, 700))
    ax1 = Axis(fig1[1, 1], xlabel = "PROD bhp limit, period 1 (bar)",
        ylabel = "NPV (million USD)",
        title = "NPV vs bhp limit — adjoint tangents in red")
    lines!(ax1, bs_coarse ./ bar, f_coarse ./ MUSD, color = :steelblue)
    scatter!(ax1, bs_coarse ./ bar, f_coarse ./ MUSD, color = :steelblue, markersize = 5)
    for (b, f, g) in zip(bs_adj, f_adj, g_adj)
        Δ = 4*bar
        lines!(ax1, [(b - Δ)/bar, (b + Δ)/bar],
            [(f - g*Δ)/MUSD, (f + g*Δ)/MUSD], color = :red, linewidth = 2)
    end
    scatter!(ax1, bs_adj ./ bar, f_adj ./ MUSD, color = :red, markersize = 9, marker = :diamond)
    # control-regime: number of period-1 steps on bhp at diag points
    ax1b = Axis(fig1[2, 1], xlabel = "PROD bhp limit, period 1 (bar)",
        ylabel = "# report steps on :bhp",
        title = "Operating-control pattern (period 1 of 8 steps / total of 24)")
    ds = sort(collect(keys(diag)))
    stairs!(ax1b, ds ./ bar, [diag[b].n_bhp_p1 for b in ds], color = :darkorange, label = "period 1")
    stairs!(ax1b, ds ./ bar, [diag[b].n_bhp for b in ds], color = :purple, label = "total")
    axislegend(ax1b, position = :rt)
    save(joinpath(figdir, "1_npv_global.png"), fig1)

    # Fig 2: zoom curve
    fig2 = Figure(size = (1100, 500))
    ax2 = Axis(fig2[1, 1], xlabel = "PROD bhp limit, period 1 (bar)",
        ylabel = "NPV (million USD)",
        title = "Zoom near base value: jump vs kink structure")
    scatterlines!(ax2, bs_zoom ./ bar, f_zoom ./ MUSD, color = :steelblue, markersize = 6)
    vlines!(ax2, [120.0], color = :gray, linestyle = :dash)
    for j in jump_ix
        vlines!(ax2, [0.5*(bs_zoom[j] + bs_zoom[j+1])/bar], color = :red, linestyle = :dot)
    end
    # adjoint tangent at 120
    let b = b0, f = f_adj[i120], g = g_adj[i120], Δ = 1.0*bar
        lines!(ax2, [(b - Δ)/bar, (b + Δ)/bar], [(f - g*Δ)/MUSD, (f + g*Δ)/MUSD],
            color = :red, linewidth = 2, label = "adjoint tangent @120")
    end
    axislegend(ax2, position = :rb)
    save(joinpath(figdir, "2_npv_zoom.png"), fig2)

    # Fig 3: gradient comparison — numerical slope of coarse curve vs adjoint
    fig3 = Figure(size = (1100, 500))
    ax3 = Axis(fig3[1, 1], xlabel = "PROD bhp limit, period 1 (bar)",
        ylabel = "dNPV/d(bhp) (USD per bar)",
        title = "Numerical slope of sweep curve vs adjoint gradient")
    b_mid = 0.5 .* (bs_coarse[1:end-1] .+ bs_coarse[2:end])
    slope_num = diff(f_coarse) ./ diff(bs_coarse)
    scatterlines!(ax3, b_mid ./ bar, slope_num .* bar, color = :steelblue,
        label = "sweep slope (central, 2.5 bar)")
    scatter!(ax3, bs_adj ./ bar, g_adj .* bar, color = :red, marker = :diamond,
        markersize = 11, label = "adjoint")
    hlines!(ax3, [0.0], color = :gray, linestyle = :dash)
    axislegend(ax3, position = :rt)
    save(joinpath(figdir, "3_gradient_comparison.png"), fig3)

    println("Figures written to $figdir")
else
    println("GLMakie not available - skipping plots.")
end

println("DONE")
