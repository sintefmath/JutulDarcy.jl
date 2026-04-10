# # Discrete Fracture Model: Validation
# <tags: Immiscible, Geothermal, Validation, Discretizations>
# This example validates the *discrete fracture model* (DFM) implementation
# in JutulDarcy against a full-dimensional reference solution for both a
# **two-phase immiscible** and a **geothermal** system.
#
# ## Background
# Fractures in porous media are thin, high-permeability features that
# significantly affect fluid flow and heat transport. In this example, we
# demonstrate two common approaches for modeling fractures:
#
# 1. **Full-dimensional (FD)**: The fracture is explicitly resolved in the mesh
#    as thin cells with high permeability. This is accurate but requires very
#    fine grids to resolve the fracture aperture.
#
# 2. **Discrete Fracture Modeling (DFM)**: The fracture is represented as a
#    codimension-one entity (surface in 3D, line in 2D) embedded in the matrix
#    mesh. Flow between matrix and fracture is handled via cross-terms. This
#    approach avoids the need to resolve the fracture aperture in the mesh and
#    gives much more flexibility in handling complex fracture geometries.
#
# In this example, we set up a fractured cube with an injector–producer
# doublet and compare the two approaches — first for two-phase immiscible flow
# and then for geothermal heat transport. We demonstrate that the DFM gives
# results that closely match the full-dimensional reference in both cases.

# ## Setup
using Jutul, JutulDarcy, HYPRE, GLMakie
import Jutul.CutCellMeshes: PlaneCut, cut_mesh

Darcy, bar, kg, meter, Kelvin, year = si_units(:darcy, :bar, :kilogram, :meter, :Kelvin, :year)

# ## Domain and simulation setup
# We create a cube domain with one fracture plane in each coordinate direction,
# placed at mid-span. The matrix has low permeability, while the fractures have
# permeability given by the cubic law ``k_f = a^2/12`` where ``a`` is the
# fracture aperture.
#
# The function `setup_fractured_cube` builds **two** simulation cases from the
# same physical domain:
#
# - A **full-dimensional** (FD) Cartesian grid where thin layers of cells
#   explicitly represent the fractures. This serves as the reference solution.
# - A **discrete fracture model** (DFM) cut mesh where fractures are embedded
#   surfaces. Cross-terms handle the fracture–matrix coupling.
#
# Both setups share the same injector–producer well pair along the ``z``-axis,
# placed in opposing corners of the domain. The function accepts a
# `fluid_model` keyword (`:two_phase` or `:geothermal`) to select the physics
# and returns ready-to-simulate `JutulCase` objects.
function setup_fractured_cube(n, L;
        matrix_permeability = 10e-3Darcy,
        matrix_porosity = 0.05,
        fracture_aperture = 1e-4,
        fracture_porosity = 0.5,
        intersection_strategy = :keep,
        well_direction = :z,
        pv_frac = 0.5,
        total_time = year,
        fluid_model = :two_phase,
        temperature = convert_to_si(90.0, :Celsius),
        injection_temperature = convert_to_si(10.0, :Celsius))

    ## Normalize input: allow scalar or 3-tuple/vector for n and L.
    dims0 = n isa Integer ? (n, n, n) : Tuple(n)
    Lxyz  = L isa Number  ? (L, L, L) : Tuple(L)
    @assert length(dims0) == 3 "n must be an Int or a length-3 tuple/vector."
    @assert length(Lxyz)  == 3 "L must be a Number or a length-3 tuple/vector."

    ## ── Full-dimensional (FD) mesh ──────────────────────────────────────
    ## Refine the Cartesian grid around each fracture plane so that a layer of
    ## cells with thickness equal to the aperture represents the fracture.
    fd_axis_widths = ntuple(3) do d
        nd = dims0[d]; Ld = Lxyz[d]
        xc = collect(range(0.0, Ld; length = nd + 1))
        xf = Ld / 2
        x = sort(vcat(xc, xf .- fracture_aperture / 2, xf .+ fracture_aperture / 2))
        x = unique(filter(xi -> 0.0 <= xi <= Ld, x))
        diff(x)
    end
    fd_cart_dims = map(length, fd_axis_widths)

    fd_mesh = UnstructuredMesh(CartesianMesh(fd_cart_dims, fd_axis_widths))
    fd_ijk = reinterpret(reshape, Int,
        map(c -> cell_ijk(fd_mesh, c), 1:number_of_cells(fd_mesh)))

    fd_matrix = reservoir_domain(fd_mesh;
        permeability = matrix_permeability, porosity = matrix_porosity)
    is_frac = [any(isapprox.(cell_dims(fd_mesh, c), fracture_aperture; rtol = 1e-2))
               for c in 1:number_of_cells(fd_mesh)]
    fd_matrix[:permeability][is_frac] .= fracture_aperture^2 / 12
    fd_matrix[:porosity][is_frac] .= fracture_porosity
    fd_matrix[:fracture_cells, Cells()] = is_frac

    radius = 0.1
    fd_inj = setup_well(fd_matrix,
        findall(fd_ijk[1, :] .== 1 .&& fd_ijk[2, :] .== 1);
        name = :Injector, simple_well = false, dir = well_direction, radius = radius)
    fd_prod = setup_well(fd_matrix,
        findall(fd_ijk[1, :] .== maximum(fd_ijk[1, :]) .&& fd_ijk[2, :] .== maximum(fd_ijk[2, :]));
        name = :Producer, simple_well = false, dir = well_direction, radius = radius)
    fd_pv = sum(pore_volume(fd_matrix))

    ## ── Lower-dimensional (DFM) mesh ────────────────────────────────────
    ## Cut a coarse Cartesian grid with fracture planes; the cuts produce an
    ## embedded fracture mesh handled by cross-terms.
    ld_mesh = UnstructuredMesh(CartesianMesh(dims0, Lxyz))
    ld_ijk = reinterpret(reshape, Int,
        map(c -> cell_ijk(ld_mesh, c), 1:number_of_cells(ld_mesh)))

    cuts = PlaneCut[]
    for d in 1:3
        center = zeros(3); center[d] = Lxyz[d] / 2
        normal = zeros(3); normal[d] = 1.0
        push!(cuts, PlaneCut(center, normal))
    end

    ld_cut_mesh, info = cut_mesh(ld_mesh, cuts; extra_out = true)
    fracture_faces = findall(info[:face_index] .== 0)
    fracture_mesh = Jutul.EmbeddedMesh(ld_cut_mesh, fracture_faces;
        intersection_strategy = intersection_strategy)
    ld_matrix = reservoir_domain(ld_cut_mesh;
        permeability = matrix_permeability, porosity = matrix_porosity)
    ld_fractures = JutulDarcy.fracture_domain(fracture_mesh, ld_matrix;
        aperture = fracture_aperture, porosity = fracture_porosity)
    ijk = ld_ijk[:, info[:cell_index]]

    ## Sort well cells by distance to the corner for consistent ordering.
    function sorted_well_cells(domain, ijk, d1_val, d2_val)
        cells = findall(ijk[1, :] .== d1_val .&& ijk[2, :] .== d2_val)
        x = domain[:cell_centroids][:, cells]
        cells[sortperm(vec(sum(x .^ 2; dims = 1)))]
    end

    inj_cells = sorted_well_cells(ld_matrix, ijk, 1, 1)
    ld_inj = setup_well(ld_matrix, inj_cells;
        name = :Injector, dir = well_direction, radius = radius, simple_well = false)
    prod_cells = sorted_well_cells(ld_matrix, ijk, maximum(ijk[1, :]), maximum(ijk[2, :]))
    ld_prod = setup_well(ld_matrix, prod_cells;
        name = :Producer, dir = well_direction, radius = radius, simple_well = false)

    ld_pv = sum(pore_volume(ld_matrix))
    if ld_fractures !== nothing
        ld_pv += sum(pore_volume(ld_fractures))
    end

    ## ── Fluid model selection ───────────────────────────────────────────
    thermal = false

    if fluid_model == :two_phase
        rho = (800.0kg / meter^3, 1000.0kg / meter^3)
        system = ImmiscibleSystem((AqueousPhase(), LiquidPhase());
            reference_densities = rho)
        init = Dict(:Pressure => 50bar, :Saturations => [0.0, 1.0])
        mix = [1.0, 0.0]

    elseif fluid_model == :geothermal
        rho = (1000.0kg / meter^3,)
        system = :geothermal
        init = Dict(:Pressure => 50bar, :Temperature => temperature)
        mix = [1.0]
        thermal = true

    else
        error("Unknown fluid_model = $fluid_model. ",
              "Supported: :two_phase, :geothermal")
    end

    ## ── Helper: build a JutulCase from a domain + wells + fractures ─────
    function make_case(matrix, wells, fractures, pv)
        nc = number_of_cells(physical_representation(matrix))

        if isnothing(fractures)
            model = setup_reservoir_model(matrix, system;
                thermal = thermal, wells = wells)
        else
            model = JutulDarcy.setup_fractured_reservoir_model(
                matrix, fractures, system;
                wells = wells, thermal = thermal)
        end

        nstep = 100
        dt = fill(total_time / nstep, nstep)
        rate = pv_frac * pv / sum(dt)

        inj_kwargs = Dict{Symbol, Any}(:density => rho[1])
        if thermal
            inj_kwargs[:temperature] = injection_temperature
        end
        ctrl_inj = InjectorControl(TotalRateTarget(rate), mix; inj_kwargs...)
        ctrl_prod = ProducerControl(BottomHolePressureTarget(10.0bar))

        control = Dict(:Injector => ctrl_inj, :Producer => ctrl_prod)
        forces = setup_reservoir_forces(model; control = control)

        state0 = setup_reservoir_state(model; init...)

        return JutulCase(model, dt, forces; state0 = state0)
    end

    case_fd = make_case(fd_matrix, [fd_inj, fd_prod], nothing, fd_pv)
    case_ld = make_case(ld_matrix, [ld_inj, ld_prod], ld_fractures, ld_pv)

    return (case_fd, case_ld)
end;

# ## Mapping helpers
# Map each DFM cell (matrix and fracture) to the nearest FD cell so that we
# can compute pointwise errors between the two solutions.
function map_dfm_to_fd(fd_domain, dfm_domain, dfm_fractures)
    x_fd = fd_domain[:cell_centroids]
    nearest(x) = last(findmin(vec(sum((x .- x_fd) .^ 2; dims = 1))))
    ld2fd_m = [nearest(c) for c in eachcol(dfm_domain[:cell_centroids])]
    ld2fd_f = [nearest(c) for c in eachcol(dfm_fractures[:cell_centroids])]
    return ld2fd_m, ld2fd_f
end

function reconstruct_full_state(state_m, state_f, state_fd, ld2fd_m, ld2fd_f)
    full = deepcopy(state_fd)
    nc = length(state_m[:Pressure])
    for k in keys(state_m)
        v, vm, vf = full[k], state_m[k], state_f[k]
        if vm isa Vector && length(vm) == nc
            v[ld2fd_m] .= vm
            v[ld2fd_f] .= vf
        elseif vm isa AbstractArray && size(vm, 2) == nc
            v[:, ld2fd_m] .= vm
            v[:, ld2fd_f] .= vf
        end
        full[k] = v
    end
    return full
end;

# ## Part 1 – Two-phase immiscible validation
# We first validate the DFM for a two-phase immiscible system. A liquid
# phase displaces the aqueous phase through the fractured cube. The injection
# rate is chosen so that 5 pore volumes are injected in one year.
cart_dims = (21, 21, 5)
phys_dims = (100.0, 100.0, 10.0)

case_fd_tp, case_dfm_tp = setup_fractured_cube(
    cart_dims, phys_dims;
    fluid_model = :two_phase,
    matrix_permeability = 1e-2Darcy,
    matrix_porosity = 0.01,
    fracture_aperture = 1e-4,
    pv_frac = 5.0,
    total_time = year);

# Simulate both cases.
sim_args = (info_level = 0, max_nonlinear_iterations = 15, max_timestep_cuts = 25)
res_fd_tp = simulate_reservoir(case_fd_tp; sim_args...)
res_dfm_tp = simulate_reservoir(case_dfm_tp; sim_args...);

# ### Two-phase well comparison
# If the DFM correctly reproduces the fracture physics, injection and
# production curves should closely follow the full-dimensional reference.
plot_well_results([res_fd_tp.wells, res_dfm_tp.wells];
    names = ["Full-dimensional", "DFM"])

# ### Two-phase cell-wise error
# Map the DFM solution back onto the FD grid and compute cell-wise saturation
# differences.
fd_domain_tp = case_fd_tp.model[:Reservoir].data_domain
dfm_domain_tp = case_dfm_tp.model[:Reservoir].data_domain
dfm_fractures_tp = case_dfm_tp.model[:Fractures].data_domain
ld2fd_m_tp, ld2fd_f_tp = map_dfm_to_fd(fd_domain_tp, dfm_domain_tp, dfm_fractures_tp)

states_fd_tp = map(s -> s[:Reservoir], res_fd_tp.result.states)
states_dfm_tp = map(s -> s[:Reservoir], res_dfm_tp.result.states)
states_frac_tp = map(s -> s[:Fractures], res_dfm_tp.result.states)

states_recon_tp = [
    reconstruct_full_state(sm, sf, sfd, ld2fd_m_tp, ld2fd_f_tp)
    for (sm, sf, sfd) in zip(states_dfm_tp, states_frac_tp, states_fd_tp)
]
delta_tp = JutulDarcy.delta_state(states_recon_tp, states_fd_tp)
plot_reservoir(case_fd_tp.model[:Reservoir], delta_tp; key = :Saturations)

# ## Part 2 – Geothermal validation
# Next we validate the DFM for a single-phase geothermal system. Cold water (10
# °C) is injected into a hot reservoir (90 °C), with fractures acting as
# preferential pathways for heat transport. We inject 50 pore volumes over one
# year so that the cold water has time to cool down the surrounding rock.
case_fd_gt, case_dfm_gt = setup_fractured_cube(
    cart_dims, phys_dims;
    fluid_model = :geothermal,
    matrix_permeability = 1e-3Darcy,
    matrix_porosity = 0.01,
    fracture_aperture = 1e-4,
    pv_frac = 50.0,
    total_time = year,
    temperature = convert_to_si(90.0, :Celsius),
    injection_temperature = convert_to_si(10.0, :Celsius));

# Simulate both cases.
sim_args_gt = (info_level = 0, initial_dt = 5.0)
res_fd_gt = simulate_reservoir(case_fd_gt; sim_args_gt...)
res_dfm_gt = simulate_reservoir(case_dfm_gt; sim_args_gt...);

# ### Geothermal well comparison
# The key metric is the **thermal breakthrough time**: when the cold injected
# water reaches the producer. If the DFM correctly captures the fracture flow
# paths, the production temperature should closely track the full-dimensional
# reference.
plot_well_results([res_fd_gt.wells, res_dfm_gt.wells],
    names = ["Full-dimensional", "DFM"])

# ### Geothermal cell-wise temperature error
# Map the DFM solution onto the FD grid and plot the cell-wise temperature
# difference. Small differences confirm that the DFM faithfully captures
# thermal transport through fractures.
fd_domain_gt = case_fd_gt.model[:Reservoir].data_domain
dfm_domain_gt = case_dfm_gt.model[:Reservoir].data_domain
dfm_fractures_gt = case_dfm_gt.model[:Fractures].data_domain
ld2fd_m_gt, ld2fd_f_gt = map_dfm_to_fd(fd_domain_gt, dfm_domain_gt, dfm_fractures_gt)

states_fd_gt = map(s -> s[:Reservoir], res_fd_gt.result.states)
states_dfm_gt = map(s -> s[:Reservoir], res_dfm_gt.result.states)
states_frac_gt = map(s -> s[:Fractures], res_dfm_gt.result.states)

states_recon_gt = [
    reconstruct_full_state(sm, sf, sfd, ld2fd_m_gt, ld2fd_f_gt)
    for (sm, sf, sfd) in zip(states_dfm_gt, states_frac_gt, states_fd_gt)
]
delta_gt = JutulDarcy.delta_state(states_recon_gt, states_fd_gt)
plot_reservoir(case_fd_gt.model[:Reservoir], delta_gt; key = :Temperature)
