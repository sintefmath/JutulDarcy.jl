using GLMakie
using JutulDarcy, Jutul
import JutulDarcy: table_to_relperm, PhaseRelativePermeability, brooks_corey_relperm, let_relperm, ReservoirRelativePermeabilities


function setup_model_with_relperm(kr; mu_ratio = 1.0)
    mesh = CartesianMesh(1000, 1000.0)
    domain = reservoir_domain(mesh)
    nc = number_of_cells(domain)
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()), reference_densities = [1000.0, 700.0])
    model,  = setup_reservoir_model(domain, sys)
    rmodel = reservoir_model(model)
    replace_variables!(rmodel, RelativePermeabilities = kr)
    JutulDarcy.add_relperm_parameters!(rmodel)
    parameters = setup_parameters(model)
    parameters[:Reservoir][:PhaseViscosities] = repeat([mu_ratio*1e-3, 1e-3], 1, nc)
    return (model, parameters)
end

function simulate_bl(model, parameters; pvi = 0.75)
    domain = reservoir_domain(model)
    nc = number_of_cells(domain)
    nstep = nc
    timesteps = fill(3600*24/nstep, nstep)
    p0 = 100*si_unit(:bar)
    tot_time = sum(timesteps)
    pv = pore_volume(domain)
    irate = pvi*sum(pv)/tot_time
    src  = SourceTerm(1, irate, fractional_flow = [1.0, 0.0])
    bc = FlowBoundaryCondition(nc, p0/2)
    forces = setup_reservoir_forces(model, sources = src, bc = bc)
    state0 = setup_reservoir_state(model, Pressure = p0, Saturations = [0.0, 1.0])
    ws, states = simulate_reservoir(state0, model, timesteps,
        forces = forces, parameters = parameters, info_level = -1)
    return states[end][:Saturations][1, :]
end

function define_relperm_single_region(krw_t, krow_t, sw_t)
    krow = PhaseRelativePermeability(reverse(1.0 .- sw), reverse(krow_t), label = :ow)
    krw = PhaseRelativePermeability(sw, krw_t, label = :ow)
    kr_def = ReservoirRelativePermeabilities(w = krw, ow = krow)
end

function simulate_and_plot(sw, krw, krow, name)
    fig = Figure(size = (1200, 800))
    ax = Axis(fig[1, 1], title = name)
    lines!(ax, sw, krw, label = L"K_{rw}")
    lines!(ax, sw, krow, label = L"K_{row}")
    axislegend(position = :ct)

    krdef = define_relperm_single_region(krw, krow, sw)
    ax = Axis(fig[2, 1], title = L"\text{Injection front at 0.75 PVI}")

    for mu_ratio in [10, 5, 1, 0.2, 0.1]
        println("Simulating with mu_ratio = $mu_ratio")
        model, parameters = setup_model_with_relperm(krdef, mu_ratio = mu_ratio)
        water_front = simulate_bl(model, parameters)
        lines!(ax, water_front, label = L"\mu_w/\mu_o=%$mu_ratio")
    end
    axislegend()
    fig
end

sw = range(0, 1, 100)
so = 1.0 .- sw

swi = 0.1
srow = 0.15
r_tot = swi + srow

krw = brooks_corey_relperm.(sw, n = 2, residual = swi, residual_total = r_tot)
krow = brooks_corey_relperm.(so, n = 3, residual = srow, residual_total = r_tot)

simulate_and_plot(sw, krw, krow, L"\text{Brooks-Corey} (N_w = 2, N_{ow} = 3)")
##
krw = brooks_corey_relperm.(sw, n = 4, residual = swi, residual_total = r_tot)
krow = brooks_corey_relperm.(so, n = 1.5, residual = srow, residual_total = r_tot)

simulate_and_plot(sw, krw, krow, L"\text{Brooks-Corey} (N_w = 4, N_{ow} = 1.5)")
# ## LET table as fully linear
krw = let_relperm.(sw, L = 1.0, E = 1.0, T = 1.0, residual = swi, residual_total = r_tot)
krow = let_relperm.(so, L = 1.0, E = 1.0, T = 1.0, residual = srow, residual_total = r_tot)

simulate_and_plot(sw, krw, krow, L"\text{LET relperm} (linear)")

# ## LET table, nonlinear exponents
krw = let_relperm.(sw, L = 3.0, E = 2.0, T = 1.0, residual = swi, residual_total = r_tot)
krow = let_relperm.(so, L = 2.0, E = 1.0, T = 2.0, residual = srow, residual_total = r_tot)

simulate_and_plot(sw, krw, krow, L"\text{LET relperm}")
##
using Optim

s_to_match  = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
kr_to_match = [0.0, 0.0, 0.01, 0.025, 0.1, 0.2, 0.4, 0.6, 0.90, 0.90, 0.90]

function relperm_mismatch(x)
    L, E, T, r, kr_max = x
    kr = let_relperm.(s_to_match, L = L, E = E, T = T, residual = r, kr_max = kr_max)
    return sum((kr .- kr_to_match).^2)
end

x0 = [1.0, 1.0, 1.0, 0.0, 1.0]
lower = [0.0, 0.0, 0.0, 0.0, 0.1]
upper = [20, 20, 20, 0.5, 1.0]
res = optimize(relperm_mismatch, lower, upper, x0, NelderMead(), Optim.Options(iterations = 100000))
L, E, T, r, kr_max = Optim.minimizer(res)

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], title = "Optimized LET Relative Permeability")

s_fine = range(0, 1, 100)
initial_kr = let_relperm.(s_fine, L = x0[1], E = x0[2], T = x0[3], residual = x0[4], kr_max = x0[5])
optimized_kr = let_relperm.(s_fine, L = L, E = E, T = T, residual = r, kr_max = kr_max)
scatter!(ax, s_to_match, kr_to_match, label = "Data")
lines!(ax, s_fine, optimized_kr, label = "Optimized LET")
lines!(ax, s_fine, initial_kr, label = "Initial LET")
axislegend(position = :ct)
fig

##
function simulate_and_plot_two_regions(sw, krw1, krow1, krw2, krow2, name)
    fig = Figure(size = (1200, 800))
    ax = Axis(fig[1, 1], title = name)
    lines!(ax, sw, krw, label = L"K_{rw}")
    lines!(ax, sw, krow, label = L"K_{row}")
    axislegend(position = :ct)

    krdef = define_relperm_single_region(krw, krow, sw)
    ax = Axis(fig[2, 1], title = L"\text{Injection front at 0.75 PVI}")

    for mu_ratio in [10, 5, 1, 0.2, 0.1]
        println("Simulating with mu_ratio = $mu_ratio")
        model, parameters = setup_model_with_relperm(krdef, mu_ratio = mu_ratio)
        water_front = simulate_bl(model, parameters)
        lines!(ax, water_front, label = L"\mu_w/\mu_o=%$mu_ratio")
    end
    axislegend()
    fig
end



# TODOs:
# Example with regions
# Example with SWOF table
##
LET = JutulDarcy.ParametricLETRelativePermeabilities()

sfront = simulate_bl(LET)
lines(sfront)