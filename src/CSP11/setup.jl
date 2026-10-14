function setup_case_a(dims_or_domain = (280, 1, 120); kwarg...)
    return first(setup_spe11_case(dims_or_domain; case = :a, kwarg...))
end

function setup_case_b(dims_or_domain = (840, 1, 120); kwarg...)
    return first(setup_spe11_case(dims_or_domain; case = :b, kwarg...))
end

"""
    case = CSP11.setup_case_c()
    case = CSP11.setup_case_c((168, 100, 120))

Set up SPE11C case (3D extension of case B).
"""
function setup_case_c(dims_or_domain = (168, 100, 120); kwarg...)
    return first(setup_spe11_case(dims_or_domain; case = :c, kwarg...))
end

"""
    setup_spe11_case(domain_or_dims; case = :b, ...)

Set up SPE11A, SPE11B or SPE11C from either an existing reservoir `DataDomain`
or a tuple of Cartesian dimensions. For A/B, dimensions are `(nx, nz)` (or
`(nx, 1, nz)`); for C they are `(nx, ny, nz)`. Case A is isothermal at 20 °C.

An existing domain is augmented in-place with benchmark wells, reporting
regions, observation points and boundary conditions. It must provide
`:permeability` and `:porosity`; `:satnum` is inferred from the geometry when
missing.

"""
function setup_spe11_case(domain_or_dims;
        case = :b,
        name = nothing,
        wells = nothing,
        input_data = Dict{String, Any}("spe11_case" => String(case)),
        domain_kwargs = NamedTuple(),
        well_kwargs = NamedTuple(),
        thermal = case != :a,
        nstep_initialization = ifelse(thermal, 10, 0),
        nstep_injection1 = 50,
        nstep_injection2 = 50,
        nstep_migration = 100,
        use_reporting_steps = false,
        divide_c_wells = true,
        kgrad = :tpfa,
        kwarg...
    )
    check_spe11_case(case)
    case == :a && thermal && throw(ArgumentError("SPE11A is isothermal; use thermal = false"))
    if domain_or_dims isa Tuple && all(x -> x isa Integer, domain_or_dims)
        dims = Tuple(Int.(domain_or_dims))
        domain = setup_spe11_domain(dims; case = case, domain_kwargs...)
        domain_name = "cartesian_$(join(dims, 'x'))"
    elseif domain_or_dims isa DataDomain
        domain = domain_or_dims
        prepare_spe11_domain!(domain, case)
        domain_name = "domain"
    else
        throw(ArgumentError("Expected a reservoir DataDomain or a tuple of Cartesian dimensions"))
    end

    if isnothing(wells)
        wells = setup_spe11_wells!(domain, case; divide_c_wells = divide_c_wells, well_kwargs...)
    elseif !haskey(domain, :well_cells)
        throw(ArgumentError("A domain used with custom wells must define domain[:well_cells]"))
    end
    if isnothing(name)
        name = domain_name
    end
    spec = spe11_case_specification(case)
    if !haskey(domain, :well_rates) || !haskey(domain, :num_well_cells)
        length(wells) == 2 || throw(ArgumentError("Custom wells must define :well_rates and :num_well_cells on the domain"))
        domain[:well_rates, nothing] = fill(spec.injection_rate, 2)
        domain[:num_well_cells, nothing] = [1, 1]
    end
    model, parameters = setup_reservoir_model_csp11(domain;
        case = case, wells = wells, thermal = thermal, dT_max_abs = 30.0,
        kgrad = kgrad, kwarg...)
    check_spe11_volumes(model, parameters, case)
    rates2 = copy(domain[:well_rates])
    rates1 = copy(rates2)
    rates1[domain[:num_well_cells][1]+1:end] .= 0.0
    forces, dt = setup_reservoir_forces_and_timesteps_csp11(model, case;
        well_labels = JutulDarcy.well_symbols(model),
        rate_injection1 = rates1, rate_injection2 = rates2,
        nstep_initialization = nstep_initialization,
        nstep_injection1 = nstep_injection1,
        nstep_injection2 = nstep_injection2,
        nstep_migration = nstep_migration,
        use_reporting_steps = use_reporting_steps)
    othername = thermal ? "thermal" : "isothermal"
    name = "spe11$(case)_$(name)_$(othername)_$kgrad"
    state0 = setup_state0_csp11(model, case)

    simulation_case = JutulCase(model, dt, forces;
        state0 = state0,
        parameters = parameters,
        input_data = input_data
    )
    return (simulation_case, name)
end

"""
    setup_spe11_domain(dims; case = :b, kwarg...)

Construct a complete Cartesian SPE11 reservoir domain directly from the
published geometry. `dims` is `(nx, nz)` for A/B and `(nx, ny, nz)` for C. The
case-C mesh is warped into physical space according to the benchmark mapping.
"""
function setup_spe11_domain(dims::Tuple; case = :b, kwarg...)
    check_spe11_case(case)
    all(x -> x isa Integer, dims) || throw(ArgumentError("Cartesian dimensions must be integers"))
    dims = Tuple(Int.(dims))
    all(>(0), dims) || throw(ArgumentError("All Cartesian dimensions must be positive"))
    spec = spe11_case_specification(case)
    if case in (:a, :b)
        if length(dims) == 2
            nx, nz = dims
            mesh_dims = (nx, 1, nz)
        elseif length(dims) == 3 && dims[2] == 1
            mesh_dims = dims
        else
            throw(ArgumentError("SPE11$(uppercase(string(case))) dimensions must be (nx, nz) or (nx, 1, nz)"))
        end
    else
        length(dims) == 3 || throw(ArgumentError("SPE11C dimensions must be (nx, ny, nz)"))
        mesh_dims = dims
    end

    mesh = UnstructuredMesh(CartesianMesh(mesh_dims, spec.physical_size), z_is_depth = true)
    if case == :c
        for (i, p) in enumerate(mesh.node_points)
            mesh.node_points[i] = SVector(p[1], p[2], p[3] - spe11c_elevation_offset(p[2]))
        end
    end

    cc = tpfv_geometry(mesh).cell_centroids
    satnum = Vector{Int}(undef, size(cc, 2))
    for i in eachindex(satnum)
        satnum[i] = spe11_facies(cc[1, i], cc[2, i], cc[3, i]; case = case)
    end
    permeability, porosity = rock_props_from_satnum(satnum, case)
    if case == :c
        permeability = transform_spe11c_permeability(permeability, cc[2, :])
    end
    domain = reservoir_domain_csp11(mesh, case;
        satnum = satnum,
        permeability = permeability,
        porosity = porosity,
        kwarg...
    )
    prepare_spe11_domain!(domain, case)
    return domain
end

function transform_spe11c_permeability(permeability, y)
    nc = length(y)
    size(permeability) == (3, nc) || throw(ArgumentError("Expected a 3×$nc diagonal permeability array"))
    out = zeros(eltype(permeability), 6, nc)
    for i in 1:nc
        k_x = permeability[1, i]
        k_y = permeability[2, i]
        k_v = permeability[3, i]
        # Derivative of the physical elevation offset in Eq. (4.1). Our third
        # coordinate is depth-positive, hence the negative yz cross term.
        slope = -3/25*((y[i] - 2500.0)/2500.0) + 1/500
        out[:, i] .= (k_x, 0.0, 0.0, k_y, -slope*k_y, k_v + slope^2*k_y)
    end
    return out
end

function spe11_thermal_properties(satnum, case)
    spec = spe11_case_specification(case)
    nc = length(satnum)
    diffusion = repeat(collect(spec.diffusion), 1, nc)
    diffusion[:, satnum .== 7] .= 0.0
    return (
        rock_thermal_conductivity = [spec.rock_conductivity[i] for i in satnum],
        diffusion = diffusion,
        fluid_thermal_conductivity = repeat([0.6, 0.088], 1, nc),
        rock_density = fill(spec.rock_density, nc),
        rock_heat_capacity = fill(spec.rock_heat_capacity, nc)
    )
end

function add_spe11_thermal_properties!(domain, case)
    for (key, value) in pairs(spe11_thermal_properties(domain[:satnum], case))
        if !haskey(domain, key)
            domain[key, Cells()] = value
        end
    end
    if !haskey(domain, :temperature)
        domain[:temperature, Cells()] = spe11_temperature.(domain[:cell_centroids][3, :], case)
    end
    return domain
end

function spe11_reporting_regions!(domain, case)
    haskey(domain, :A) && haskey(domain, :B) && haskey(domain, :C) && return domain
    cc = domain[:cell_centroids]
    nc = size(cc, 2)
    A = zeros(nc)
    B = zeros(nc)
    C = zeros(nc)
    boxes = spe11_case_specification(case).reporting_boxes
    for i in 1:nc
        x = cc[1, i]
        elevation = (case == :a ? 1.2 : spe11_size_z) - cc[3, i]
        for (weights, box) in zip((A, B, C), boxes)
            xmin, xmax, zmin, zmax = box
            weights[i] = xmin <= x <= xmax && zmin <= elevation <= zmax
        end
    end
    domain[:A, Cells()] = A
    domain[:B, Cells()] = B
    domain[:C, Cells()] = C
    return domain
end

function spe11_buffer_volumes(domain, case)
    check_spe11_case(case)
    satnum = domain[:satnum]
    buffer = zeros(eltype(domain[:volumes]), length(satnum))
    case == :a && return buffer
    boundary_centroids = domain[:boundary_centroids]
    isempty(domain[:boundary_neighbors]) && return buffer
    xmin, xmax = extrema(boundary_centroids[1, :])
    ymin, ymax = extrema(boundary_centroids[2, :])
    xtol = max(1.0, xmax - xmin)*1e-8
    ytol = max(1.0, ymax - ymin)*1e-8
    facies = spe11_case_specification(case).buffer_facies
    for face in eachindex(domain[:boundary_neighbors])
        cell = domain[:boundary_neighbors][face]
        x, y = boundary_centroids[1, face], boundary_centroids[2, face]
        on_x_side = abs(x - xmin) <= xtol || abs(x - xmax) <= xtol
        on_c_front_or_back = case == :c && (abs(y - ymin) <= ytol || abs(y - ymax) <= ytol)
        if (on_x_side || on_c_front_or_back) && satnum[cell] in facies
            # Add every eligible face, including both faces of corner cells.
            buffer[cell] += spe11_buffer_length*domain[:boundary_areas][face]
        end
    end
    return buffer
end

function spe11_boundary_conditions!(domain, case)
    if !haskey(domain, :spe11_physical_volumes)
        domain[:spe11_physical_volumes, Cells()] = tpfv_geometry(physical_representation(domain)).volumes
    end
    buffer = spe11_buffer_volumes(domain, case)
    domain[:spe11_buffer_volumes, Cells()] = buffer
    domain[:volumes, Cells()] = domain[:spe11_physical_volumes] .+ buffer
    is_boundary = buffer .> 0
    domain[:boundary, nothing] = findall(is_boundary)
    domain[:is_boundary, Cells()] = is_boundary

    thermal_faces = Int[]
    pressure_faces = Int[]
    for face in eachindex(domain[:boundary_neighbors])
        normal = domain[:boundary_normals][:, face]
        if case == :a && normal[3] < -0.5
            push!(pressure_faces, face)
        elseif case != :a && abs(normal[3]) > abs(normal[1]) + abs(normal[2])
            push!(thermal_faces, face)
        end
    end
    domain[:spe11_thermal_boundary_faces, nothing] = thermal_faces
    domain[:spe11_pressure_boundary_faces, nothing] = pressure_faces
    return domain
end

function spe11_pressure_boundary_conditions(model)
    domain = reservoir_domain(model)
    trans = compute_boundary_trans(domain)
    table = spe11_hydrostatic_pressure_table(model, :a)
    rho_water = reservoir_model(model).secondary_variables[:PhaseMassDensities].tab(1.1e5, 293.15)[1]
    return map(domain[:spe11_pressure_boundary_faces]) do face
        cell = domain[:boundary_neighbors][face]
        # FlowBoundaryCondition compares pressures at the cell centre. Lift the
        # prescribed face pressure hydrostatically to avoid draining a resting
        # water column through the top row.
        pressure = table(domain[:cell_centroids][3, cell])
        FlowBoundaryCondition(cell, pressure, 293.15;
            trans_flow = trans[face], fractional_flow = [1.0, 0.0], density = rho_water)
    end
end

function spe11_thermal_boundary_conditions(model; well_rates = Dict{Symbol, Float64}(),
        injection_temperature = 283.15)
    rmodel = JutulDarcy.reservoir_model(model)
    JutulDarcy.model_is_thermal(rmodel) || return nothing
    domain = reservoir_domain(model)
    phi = domain[:porosity]
    # The impermeable top/bottom facies contain water. Match the reservoir's
    # parallel rock/fluid conductivity, using the physical boundary geometry.
    fluid = domain[:fluid_thermal_conductivity]
    water = fluid isa AbstractVector ? fluid : vec(fluid[1, :])
    conductivity = (1.0 .- phi).*domain[:rock_thermal_conductivity] .+ phi.*water
    trans = compute_boundary_trans(domain, conductivity)
    # Aggregate by cell so even a one-layer grid has one force per cell. This
    # avoids concurrent writes when forces are evaluated by device kernels.
    conductance = Dict{Int, Float64}()
    weighted_temperature = Dict{Int, Float64}()
    function add_boundary!(cell, t, temperature)
        conductance[cell] = get(conductance, cell, 0.0) + t
        weighted_temperature[cell] = get(weighted_temperature, cell, 0.0) + t*temperature
    end
    for face in domain[:spe11_thermal_boundary_faces]
        cell = domain[:boundary_neighbors][face]
        t = trans[face]
        temperature = spe11_geothermal_temperature(domain[:boundary_centroids][3, face])
        add_boundary!(cell, t, temperature)
    end
    for (name, connections) in get(domain, :spe11_well_thermal_boundaries, Dict())
        if get(well_rates, name, 0.0) > 0.0
            for (cell, t) in connections
                add_boundary!(cell, t, injection_temperature)
            end
        end
    end
    return [FlowBoundaryCondition(cell, 3e7, weighted_temperature[cell]/conductance[cell];
        trans_flow = 0.0, trans_thermal = conductance[cell])
        for cell in sort!(collect(keys(conductance))) if conductance[cell] > 0.0]
end

"""
    check_spe11_volumes(model, parameters, case)

Check physical geometry, lateral buffer storage and initialized model pore
volume. Returns measured and integrated reference volumes in cubic metres.
Facies discretization differences are logged at debug level; inconsistent
geometry or model storage produces warnings with actual and expected values.
"""
function check_spe11_volumes(model, parameters, case)
    domain = reservoir_domain(model)
    reference = spe11_reference_volumes(case)
    physical = sum(domain[:spe11_physical_volumes])
    expected_bulk = physical + sum(spe11_buffer_volumes(domain, case))
    bulk = sum(domain[:volumes])
    expected_pore = sum(pore_volume(domain))
    pore = sum(pore_volume(model, parameters))
    function check(message, actual, expected; rtol = 1e-8)
        relative_error = abs(actual - expected)/expected
        if !isfinite(relative_error) || relative_error > rtol
            @warn message actual expected relative_error
        end
    end
    check("SPE11 physical bulk volume differs from the specified domain dimensions",
        physical, reference.physical_bulk)
    check("SPE11 effective bulk volume is inconsistent with the mesh and lateral buffer areas",
        bulk, expected_bulk)
    check("SPE11 model pore volume is inconsistent with the reservoir domain",
        pore, expected_pore)
    bulk_error = abs(bulk - reference.bulk)/reference.bulk
    pore_error = abs(pore - reference.pore)/reference.pore
    if max(bulk_error, pore_error) > 1e-3
        # Sampling facies at cell centres changes their volumes on coarse
        # meshes. Report this discretization error without a setup warning.
        @debug "SPE11 discretized storage volumes differ from the integrated facies geometry" bulk pore reference_bulk = reference.bulk reference_pore = reference.pore bulk_relative_error = bulk_error pore_relative_error = pore_error
    end
    return (physical = physical, bulk = bulk, pore = pore, reference = reference,
        bulk_relative_error = bulk_error, pore_relative_error = pore_error)
end

function spe11_observation_points!(domain, case)
    haskey(domain, :spe11_observation_cells) && return domain
    points = spe11_case_specification(case).observation_points
    cc = domain[:cell_centroids]
    cell_points = eachcol(cc)
    observations = zeros(Int, number_of_cells(domain))
    cells = map(point -> find_closest_point(cell_points, point), points)
    for (i, cell) in enumerate(cells)
        observations[cell] = i
    end
    domain[:spe11_observation_cells, nothing] = cells
    domain[:observation_points, Cells()] = observations
    return domain
end

function prepare_spe11_domain!(domain::DataDomain, case)
    check_spe11_case(case)
    get(domain, :spe11_case, case) == case || throw(ArgumentError("The domain was prepared for a different SPE11 case"))
    domain[:spe11_case, nothing] = case
    for key in (:permeability, :porosity)
        haskey(domain, key) || throw(ArgumentError("The supplied domain must define :$key"))
    end
    if !haskey(domain, :satnum)
        cc = domain[:cell_centroids]
        satnum = [spe11_facies(cc[1, i], cc[2, i], cc[3, i]; case = case) for i in axes(cc, 2)]
        domain[:satnum, Cells()] = satnum
    end
    domain[:porosity][domain[:satnum] .== 7] .= spe11_facies7_porosity
    add_spe11_thermal_properties!(domain, case)
    spe11_reporting_regions!(domain, case)
    spe11_boundary_conditions!(domain, case)
    spe11_observation_points!(domain, case)
    return domain
end

function spe11_trajectory_cells(domain, trajectory; n = 501)
    mesh = physical_representation(domain) |> UnstructuredMesh
    cells, extra = Jutul.find_enclosing_cells(mesh, trajectory; n = n, extra_out = true)
    isempty(cells) && error("The SPE11 well trajectory did not intersect the supplied domain")
    lengths = extra[:lengths]
    return cells, lengths
end

function setup_spe11_wells!(domain::DataDomain, case; divide_c_wells = true, kwarg...)
    spec = spe11_case_specification(case)
    # The benchmark prescribes fluxes at the holes, without a wellbore storage
    # volume. Keep a small positive native volume to avoid singular equations
    # while making well accumulation negligible, also at the lab-scale rates.
    default_options = (simple_well = true, radius = spec.well_radius, dir = :y,
        volume_multiplier = 1e-8)
    options = merge(default_options, values(kwarg))
    if case in (:a, :b)
        cc = domain[:cell_centroids]
        points = eachcol(cc)
        well_points = spec.well_positions
        well_cells = [find_closest_point(points, p) for p in well_points]
        wells = [setup_well(domain, well_cells[i]; options..., name = Symbol(:INJ, i-1)) for i in 1:2]
        domain[:well_cells, nothing] = well_cells
        domain[:num_well_cells, nothing] = [1, 1]
        domain[:well_rates, nothing] = fill(spec.injection_rate, 2)
    else
        trajectory1 = [2700.0 1000.0 900.0; 2700.0 4000.0 900.0]
        y = collect(range(1000.0, 4000.0, length = 101))
        trajectory2 = hcat(fill(5100.0, length(y)), y, 500.0 .- spe11c_elevation_offset.(y))
        cells1, lengths1 = spe11_trajectory_cells(domain, trajectory1; n = 1001)
        cells2, lengths2 = spe11_trajectory_cells(domain, trajectory2; n = 11)
        cells = [cells1; cells2]
        if divide_c_wells
            wells = Any[]
            for i in eachindex(cells)
                push!(wells, setup_well(domain, cells[i]; options..., name = Symbol(:INJ, i-1)))
            end
            rates1 = spec.injection_rate.*lengths1./sum(lengths1)
            rates2 = spec.injection_rate.*lengths2./sum(lengths2)
            n1 = length(cells1)
            n2 = length(cells2)
        else
            wells = [setup_well(domain, cells1; options..., name = :INJ0),
                setup_well(domain, cells2; options..., name = :INJ1)]
            rates1 = [spec.injection_rate]
            rates2 = [spec.injection_rate]
            n1 = n2 = 1
        end
        domain[:well_cells, nothing] = cells
        domain[:num_well_cells, nothing] = [n1, n2]
        domain[:well_rates, nothing] = [rates1; rates2].*si_unit(:kilogram)./si_unit(:second)
    end
    return wells
end

function reservoir_domain_csp11(G, case = :b; satnum,
        temperature = nothing,
        permeability = first(rock_props_from_satnum(satnum, case)),
        porosity = last(rock_props_from_satnum(satnum, case)), kwarg...)
    check_spe11_case(case)
    nc = number_of_cells(G)
    length(satnum) == nc || throw(ArgumentError("satnum must have number of cells entries ($nc), was $(length(satnum))"))
    all(in(1:7), satnum) || throw(ArgumentError("SPE11 saturation regions must be in 1:7"))
    porosity = porosity isa Real ? fill(porosity, nc) : copy(porosity)
    porosity[satnum .== 7] .= spe11_facies7_porosity
    options = merge(spe11_thermal_properties(satnum, case), values(kwarg))
    domain = reservoir_domain(G; satnum = satnum, permeability = permeability,
        porosity = porosity, options...)
    if isnothing(temperature)
        temperature = spe11_temperature.(domain[:cell_centroids][3, :], case)
    end
    domain[:temperature, Cells()] = temperature isa Real ? fill(temperature, nc) : temperature
    domain[:spe11_case, nothing] = case
    return domain
end

function rock_props_from_satnum(satnum, case)
    spec = spe11_case_specification(case)
    all(x -> x isa Integer && x in 1:7, satnum) || throw(ArgumentError("SPE11 saturation regions must be integers in 1:7"))
    perm = zeros(3, length(satnum))
    poro = zeros(length(satnum))
    for (i, reg) in enumerate(satnum)
        perm[1, i] = perm[2, i] = spec.permeability[reg]
        perm[3, i] = spec.vertical_permeability_ratio*spec.permeability[reg]
        poro[i] = spec.porosity[reg]
    end
    return perm, poro
end

function spe11_phase_relative_permeability(immobile, label;
        n = 65, exponent = spe11_case_specification(:b).exponent)
    effective_saturation = range(0.0, 1.0; length = n)
    saturation = immobile .+ (1.0 - immobile).*effective_saturation
    kr = effective_saturation.^exponent
    return PhaseRelativePermeability(saturation, kr; label = label)
end

function spe11_erf(x)
    ax = abs(x)
    if ax < 0.5
        # The power series is accurate and well-conditioned near zero.
        term = total = ax
        for n in 1:100
            term *= -ax^2/n
            increment = term/(2n + 1)
            total += increment
            abs(increment) <= eps(Float64)*abs(total) && break
        end
        value = 2.0/sqrt(pi)*total
    else
        # Compact approximation with a maximum absolute error around 3e-8.
        t = 1.0/(1.0 + 0.5ax)
        polynomial = 0.17087277
        polynomial = -0.82215223 + t*polynomial
        polynomial = 1.48851587 + t*polynomial
        polynomial = -1.13520398 + t*polynomial
        polynomial = 0.27886807 + t*polynomial
        polynomial = -0.18628806 + t*polynomial
        polynomial = 0.09678418 + t*polynomial
        polynomial = 0.37409196 + t*polynomial
        polynomial = 1.00002368 + t*polynomial
        tail = t*exp(-ax^2 - 1.26551223 + t*polynomial)
        value = 1.0 - tail
    end
    return signbit(x) ? -value : value
end

function spe11_capillary_pressure(sg, immobile, entry_pressure;
        exponent = spe11_case_specification(:b).exponent,
        max_pressure = spe11_case_specification(:b).max_capillary_pressure)
    sw = 1.0 - sg
    normalized_sw = max((sw - immobile)/(1.0 - immobile), 0.0)
    if iszero(normalized_sw)
        return max_pressure
    end
    pc_unbounded = entry_pressure*normalized_sw^(-1.0/exponent)
    erf_argument = pc_unbounded/max_pressure*sqrt(pi)/2.0
    return max_pressure*spe11_erf(erf_argument)
end

function spe11_capillary_pressure_table(immobile, entry_pressure;
        relative_tolerance = 1e-3,
        absolute_tolerance = 1.0,
        kwarg...
    )
    pc(sg) = spe11_capillary_pressure(sg, immobile, entry_pressure; kwarg...)
    effective_saturation = vcat(
        collect(range(0.0, 1.0; length = 17)),
        10.0.^range(-14.0, 0.0; length = 33)
    )
    sort!(effective_saturation)
    unique!(effective_saturation)
    points = 1.0 .- (immobile .+ (1.0 - immobile).*effective_saturation)
    push!(points, 1.0)
    sort!(points)
    unique!(points)

    function refine_interval!(lo, hi, pc_lo, pc_hi, level = 0)
        level == 40 && return
        midpoint = (lo + hi)/2.0
        pc_midpoint = pc(midpoint)
        needs_refinement = any((0.25, 0.5, 0.75)) do fraction
            sg = lo + fraction*(hi - lo)
            pc_exact = pc(sg)
            pc_linear = pc_lo + fraction*(pc_hi - pc_lo)
            tolerance = max(absolute_tolerance, relative_tolerance*abs(pc_exact))
            abs(pc_exact - pc_linear) > tolerance
        end
        if needs_refinement
            push!(points, midpoint)
            refine_interval!(lo, midpoint, pc_lo, pc_midpoint, level + 1)
            refine_interval!(midpoint, hi, pc_midpoint, pc_hi, level + 1)
        end
        return
    end

    initial_points = copy(points)
    for i in 1:(length(initial_points) - 1)
        lo = initial_points[i]
        hi = initial_points[i + 1]
        refine_interval!(lo, hi, pc(lo), pc(hi))
    end
    sort!(points)
    pressure = pc.(points)
    return Jutul.LinearInterpolant(points, pressure; constant_dx = false)
end

"""
    spe11_saturation_functions(satnum; case = :b)

Construct the SPE11 Brooks-Corey relative-permeability and capillary-
pressure functions directly from the benchmark parameters. The original
PYOPMSPE11 deck represented the same functions with 200,000 coupled `SGOF`
rows per permeable facies. Relative permeability is sampled independently on a
small uniform effective-saturation grid, while capillary pressure uses adaptive
sampling around its much sharper transition.
"""
function spe11_saturation_functions(satnum; case = :b)
    check_spe11_case(case)
    spec = spe11_case_specification(case)
    exponent = spec.exponent
    max_pressure = spec.max_capillary_pressure
    regions = Int.(vec(satnum))
    all(in(1:7), regions) || throw(ArgumentError("SPE11 saturation regions must be in 1:7"))

    gas = spe11_phase_relative_permeability(
        spe11_nonwetting_immobile_saturation, :g; exponent = exponent)
    impermeable_gas = PhaseRelativePermeability([0.0, 1.0], [0.0, 1.0]; label = :g)
    gas_tables = ntuple(i -> i == 7 ? impermeable_gas : gas, 7)

    impermeable_liquid = PhaseRelativePermeability([0.0, 1.0], [0.0, 1.0]; label = :og)
    liquid_tables = ntuple(7) do i
        if i == 7
            impermeable_liquid
        else
            spe11_phase_relative_permeability(spe11_wetting_immobile_saturation[i], :og;
                exponent = exponent)
        end
    end
    relative_permeabilities = JutulDarcy.ReservoirRelativePermeabilities(;
        og = liquid_tables,
        g = gas_tables,
        regions = regions
    )

    capillary_tables = ntuple(7) do i
        if i == 7
            # The deck used a two-point fallback for the impermeable facies.
            Jutul.LinearInterpolant(
                [0.0, 1.0], [0.0, max_pressure];
                constant_dx = false
            )
        else
            entry_pressure = spec.entry_pressure[i]
            spe11_capillary_pressure_table(
                spe11_wetting_immobile_saturation[i], entry_pressure;
                exponent = exponent, max_pressure = max_pressure)
        end
    end
    capillary_pressure = JutulDarcy.SimpleCapillaryPressure(
        (capillary_tables, ); regions = regions)
    return relative_permeabilities, capillary_pressure
end

function setup_reservoir_model_csp11(reservoir::DataDomain;
        case = get(reservoir, :spe11_case, :b), include_satfun = true,
        thermal = false, kwarg...)
    check_spe11_case(case)
    case == :a && thermal && throw(ArgumentError("SPE11A is isothermal; use thermal = false"))
    reservoir[:gdz, Faces()] = JutulDarcy.compute_face_gdz(reservoir[:neighbors], reservoir[:cell_centroids][3, :]; g = spe11_gravity)
    model, parameters = setup_reservoir_model(reservoir, :co2brine;
        co2_source = :csp11,
        thermal = thermal,
        extra_out = true,
        kwarg...
    )
    for name in JutulDarcy.well_symbols(model)
        parameters[name][:PerforationGravityDifference] .*= spe11_gravity/Jutul.gravity_constant
    end
    if include_satfun
        kr, pc = spe11_saturation_functions(reservoir[:satnum]; case = case)
        set_secondary_variables!(model[:Reservoir],
            RelativePermeabilities = kr,
            CapillaryPressure = pc
        )
    end
    if thermal
        setup_spe11_well_thermal_boundaries!(model, parameters)
    end
    return (model, parameters)
end

function setup_spe11_well_thermal_boundaries!(model, parameters)
    connections = Dict{Symbol, Vector{Tuple{Int, Float64}}}()
    for name in JutulDarcy.well_symbols(model)
        well = physical_representation(model[name].data_domain)
        conductance = parameters[name][:WellIndicesThermal]
        connections[name] = collect(zip(well.perforations.reservoir, Float64.(conductance)))
        # Apply prescribed temperatures with native reservoir thermal boundaries,
        # so disabled wells cannot conduct heat between their perforations.
        conductance .= 0.0
    end
    reservoir_domain(model)[:spe11_well_thermal_boundaries, nothing] = connections
    return model
end

function setup_reservoir_forces_and_timesteps_csp11(model, case = :b;
        well_labels = JutulDarcy.well_symbols(model),
        nstep_initialization = case == :a ? 0 : 1,
        nstep_injection1 = 1, nstep_injection2 = 1, nstep_migration = 1,
        rate_injection1 = [i == 1 ? spe11_case_specification(case).injection_rate : 0.0 for i in eachindex(well_labels)],
        rate_injection2 = fill(spe11_case_specification(case).injection_rate, length(well_labels)),
        use_reporting_steps = false,
        injection_temperature = spe11_case_specification(case).injection_temperature)
    spec = spe11_case_specification(case)
    counts = (nstep_initialization, nstep_injection1, nstep_injection2, nstep_migration)
    all(n -> n isa Integer && n >= 0, counts) || throw(ArgumentError("Step counts must be non-negative integers"))
    all(>(0), counts[2:end]) || throw(ArgumentError("Injection and migration periods must have at least one step"))
    length(well_labels) == length(rate_injection1) == length(rate_injection2) ||
        throw(ArgumentError("Each injection period must specify one rate per well"))
    all(r -> isfinite(r) && r >= 0.0, Iterators.flatten((rate_injection1, rate_injection2))) ||
        throw(ArgumentError("Injection rates must be finite and non-negative"))
    _, rho_co2 = JutulDarcy.reference_densities(reservoir_model(model).system)
    dt = Float64[]
    forces = Dict{Symbol, Any}[]
    base_bc = case == :a ? spe11_pressure_boundary_conditions(model) : spe11_thermal_boundary_conditions(model)
    function forcing(rates)
        controls = Dict(name => rate_to_injection_control(rate, rho_co2, injection_temperature)
            for (name, rate) in zip(well_labels, rates))
        bc = case == :a ? base_bc : spe11_thermal_boundary_conditions(model;
            well_rates = Dict(zip(well_labels, rates)), injection_temperature = injection_temperature)
        well_forces = Dict(name => setup_forces(model[name];
            mask = PerforationMask(fill(rate > 0.0 ? 1.0 : 0.0,
                length(physical_representation(model[name].data_domain).perforations.reservoir))))
            for (name, rate) in zip(well_labels, rates))
        return setup_reservoir_forces(model; control = controls, bc = bc, well_forces...)
    end
    disabled = forcing(zeros(length(well_labels)))
    initialization = nstep_initialization == 0 ? 0.0 : spec.initialization_time
    reservoir_domain(model)[:spe11_initialization_time, nothing] = initialization
    if initialization > 0.0
        steps = use_reporting_steps ? 1 : nstep_initialization
        add_timesteps_and_forces!(dt, forces, fill(initialization/steps, steps), disabled)
    end
    if use_reporting_steps
        times = spe11_reporting_times(case; kind = case == :a ? :sparse : :dense)
    else
        t1, t2 = spec.injection_times
        times = [0.0;
            collect(range(0.0, t1; length = nstep_injection1 + 1))[2:end];
            collect(range(t1, t1 + t2; length = nstep_injection2 + 1))[2:end];
            collect(range(t1 + t2, spec.simulation_time; length = nstep_migration + 1))[2:end]]
    end
    f1, f2 = forcing(rate_injection1), forcing(rate_injection2)
    injection1_end, injection_end = spec.injection_times[1], sum(spec.injection_times)
    for i in 2:length(times)
        start = times[i-1]
        force = start < injection1_end ? f1 : start < injection_end ? f2 : disabled
        add_timesteps_and_forces!(dt, forces, (times[i] - start,), force)
    end
    return forces, dt
end

function add_timesteps_and_forces!(dt, forces, steps, forcing)
    for step in steps
        isfinite(step) && step > 0.0 || throw(ArgumentError("Time steps must be finite and positive"))
        push!(dt, step)
        push!(forces, deepcopy(forcing))
    end
end

function rate_to_injection_control(wrate, co2_density, co2_temp; kwarg...)
    if wrate < 0
        error("Only injection supported")
    elseif wrate == 0
        ctrl = DisabledControl()
    else
        t = TotalRateTarget(wrate/co2_density)
        mix = [0.0, 1.0]
        ctrl = InjectorControl(t, mix; density = co2_density, temperature = co2_temp, kwarg...)
    end
    return ctrl
end

function spe11_hydrostatic_pressure_table(model, case)
    spec = spe11_case_specification(case)
    domain = reservoir_domain(model)
    density = reservoir_model(model).secondary_variables[:PhaseMassDensities].tab
    zmin, zmax = extrema(domain[:cell_centroids][3, :])
    zmin = min(zmin, minimum(domain[:boundary_centroids][3, :]), spec.datum_depth)
    zmax = max(zmax, maximum(domain[:boundary_centroids][3, :]), spec.datum_depth)
    density_function(p, depth, phase) = density(p, spe11_temperature(depth, case))[1]
    return JutulDarcy.phase_pressure_depth_table(spec.datum_depth, zmin, zmax,
        spec.datum_pressure, density_function, 1; g = spe11_gravity)
end

function setup_state0_csp11(model, case = get(reservoir_domain(model), :spe11_case, :b))
    domain = reservoir_domain(model)
    z = domain[:cell_centroids][3, :]
    p = spe11_hydrostatic_pressure_table(model, case).(z)
    temperature = spe11_temperature.(z, case)
    return setup_reservoir_state(model; OverallMoleFractions = [1.0, 0.0],
        Pressure = p, Temperature = temperature)
end
