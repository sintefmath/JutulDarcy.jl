# # Consistent discretizations: Average MPFA and nonlinear TPFA
# This example demonstrates how to use alternative discretizations for the
# pressure gradient term in the Darcy equation, i.e. the approximation of the
# Darcy flux:
#
# ``\mathbf{}{K}(\nabla p + \rho g \Delta z)``.
#
# It is well-known that for certain combinations of grid geometry and
# permeability fields, the classical two-point flux approximation scheme can
# give incorrect results. This is due to the fact that the TPFA scheme is not a
# formally consistent method when the product of the permeability tensor and the
# normal vector does not align with the cell-to-cell vectors over a face (lack
# of K-orthogonality).
#
# In such cases, it is often beneficial to use a consistent discretization.
# JutulDarcy includes a class of linear and nonlinear schemes that are designed
# to be accurate even for challenging grids.
# 
# For further details on this class of methods, which differ a bit from the
# classical MPFA-O type method often seen in the literature, see
# [schneider_nonlinear](@cite), [zhang_nonlinear](@cite) and
# [raynaud_discretization](@cite).

# ## Define a mesh and twist the nodes
# This makes the mesh non K-orthogonal and will lead to wrong solutions for the
# default TPFA scheme.
using Jutul
using JutulDarcy
using LinearAlgebra
using GLMakie
using Test # hide

sys = SinglePhaseSystem()
nx = nz = 100
pdims = (1.0, 1.0)
g = CartesianMesh((nx, nz), pdims)

g = UnstructuredMesh(g)
D = dim(g)

v = 0.1
for i in eachindex(g.node_points)
    x, y = g.node_points[i]
    shiftx =  v*sin(π*x)*sin(3*(-π/2 + π*y))
    shifty =  v*sin(π*y)*sin(3*(-π/2 + π*x))
    g.node_points[i] += [shiftx, shifty]
end

nc = number_of_cells(g)
domain = reservoir_domain(g, permeability = 0.1*si_unit(:darcy))

fig = Figure()
Jutul.plot_mesh_edges!(Axis(fig[1, 1]), g)
fig
# ## Create a test problem function
# We set up a problem for our given domain with left and right boundary boundary
# conditions that correspond to a linear pressure drop. We can expect the
# steady-state pressure solution to be linear between the two faces as there is
# no variation in permeability or significant compressibility. The function will
# return the pressure solution at the end of the simulation for a given scheme.
function solve_test_problem(scheme)
    model, parameters = setup_reservoir_model(domain, sys,
        general_ad = true,
        kgrad = scheme,
        block_backend = false
    )
    state0 = setup_reservoir_state(model, Pressure = 1e5)
    nc = number_of_cells(g)
    bcells = Int64[]
    bpres = Float64[]
    for k in 1:nz
        bnd_l = Jutul.cell_index(g, (1, k))
        bnd_r = Jutul.cell_index(g, (nx, k))
        push!(bcells, bnd_l)
        push!(bpres, 1e5)
        push!(bcells, bnd_r)
        push!(bpres, 2e5)
    end
    bc = flow_boundary_condition(bcells, domain, bpres)
    forces = setup_reservoir_forces(model, bc = bc)

    dt = [si_unit(:day)]
    _, states = simulate_reservoir(state0, model, dt,
        forces = forces, failure_cuts_timestep = false,
        tol_cnv = 1e-6,
        linear_solver = GenericKrylov(preconditioner = AMGPreconditioner(:smoothed_aggregation), rtol = 1e-6)
        )
    return states[end][:Pressure]
end
# ## Solve the test problem with three different schemes
# - TPFA (two-point flux approximation, inconsistent, linear)
# - Average MPFA (consistent, linear)
# - NTPFA (nonlinear two-point flux approximation, consistent, nonlinear)
results = Dict()
for m in [:tpfa, :avgmpfa, :ntpfa]
    println("Solving $m")
    results[m] = solve_test_problem(m);
end
# ## Plot the results
# We plot the pressure solution for each of the schemes, as well as the error.
# Note that the color axis varies between error plots. As the grid is quite
# skewed, we observe significant errors for the TPFA scheme, with no significant
# error for the consistent schemes.
x = domain[:cell_centroids][1, :]
get_ref(x) = x*1e5 + 1e5
x_distinct = sort(unique(x))
sol = get_ref.(x_distinct)

fig = Figure(size = (1200, 600))
for (i, (m, s)) in enumerate(results)
    ref = get_ref.(x)
    err = norm(s .- ref)/norm(ref)

    ax = Axis(fig[1, i], title = "$m")
    lines!(ax, x_distinct, sol, color = :red)
    scatter!(ax, x, s, label = m, markersize = 1, alpha = 0.5, transparency = true)

    ax2 = Axis(fig[2, i], title = "Error=$err")
    Δ = s - ref
    emin, emax = extrema(Δ)
    largest = max(abs(emin), abs(emax))
    crange = (-largest, largest)
    plt = plot_cell_data!(ax2, g, Δ, colorrange = crange, colormap = :seismic)
    Colorbar(fig[3, i], plt, vertical = false)
    if m != :tpfa  # hide
        @test err < 1e-4 # hide
    end  # hide
end
fig
