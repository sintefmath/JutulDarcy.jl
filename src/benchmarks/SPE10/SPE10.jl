module SPE10
    using LazyArtifacts, Artifacts, JutulDarcy, DelimitedFiles, Jutul

    """
        JutulDarcy.SPE10.setup_case()
        JutulDarcy.SPE10.setup_case(layers = 1:5, num_steps = 200)
        JutulDarcy.SPE10.setup_case(layers = 1:5, tiling = (2, 3, 1))

    Set up the SPE10, model2, benchmark case. This model is a 1.1 million cell
    nearly incompressible oil-water model with highly heterogenous static
    properties (permeability/porosity).

    `tiling = (tilex, tiley, tilez)` repeats the selected layers and static
    properties in each direction. Each horizontal tile has an injector at its
    center, with producers shared at tile corners to avoid adjacent wells at
    tile boundaries. Vertical wells span all repeated layers.

    Pass `system` to use a custom fluid system without installing the SPE10
    relative permeability or PVT functions. `equil` can be an `EquilibriumRegion`
    or a function `model -> EquilibriumRegion(model, ...)`. `forces_arg` forwards
    keyword arguments to `setup_well_forces`, including a custom `injection_control`.
    """
    function setup_case(;
            num_steps = 100,
            t_total = 2000*si_unit(:day),
            equil = missing,
            forces_arg = NamedTuple(),
            kwarg...
        )
        model = setup_model(; kwarg...)
        state0 = setup_state0(model; equil = equil)
        dt = fill(t_total/num_steps, num_steps)
        forces = setup_well_forces(model; forces_arg...)
        return JutulCase(model, dt, forces, state0 = state0)
    end

    function setup_reservoir(;
            layers = 1:85,
            anisotropy = true,
            minporo = 0.01,
            remove_cells = true,
            tiling = (1, 1, 1),
            source = spe10_artifact_path()
        )
        check_tiling(tiling)
        minporo > 0 || error("Minimum porosity must be positive.")
        nx = 60
        ny = 220
        nz_data = 85
        nvals_data = nx*ny*nz_data
        perm = readdlm(joinpath(source, "spe_perm.dat"))'
        poro = readdlm(joinpath(source, "spe_phi.dat"))'
        length(poro) == nvals_data || error("Porosity data must have $nvals_data values")
        length(perm) == 3*nvals_data || error("Permeability data must have 3*$nvals_data values")

        function extract_layers_and_flatten(x)
            x = vec(x)
            x = reshape(x, (nx, ny, nz_data))
            if layers != 1:85
                x = x[:, :, layers]
            end
            return vec(repeat(x, tiling...))
        end

        function get_perm(i)
            @assert size(perm) == (6, 561000)
            md = si_unit(:milli)*si_unit(:darcy)
            offset = nvals_data ÷ 6
            vals = perm[:, ((i-1)*offset + 1):(i*offset)] .* md
            vals = vec(vals)
            return extract_layers_and_flatten(vals)
        end
        poro = extract_layers_and_flatten(poro)
        permx = get_perm(1)
        permy = get_perm(2)
        permz = get_perm(3)
        if anisotropy isa Bool
            if anisotropy
                perm = cat(permx', permy', permz'; dims = 1)
            else
                perm = permx
            end
        else
            anisotropy isa Number
            permz = permx .* anisotropy
            perm = cat(permx', permx', permz'; dims = 1)
        end
        mesh = setup_mesh(layers; tiling = tiling)
        if remove_cells
            # Remove cells with too low porosity
            active = findall(poro .>= minporo)
            perm = perm isa AbstractVector ? perm[active] : perm[:, active]
            poro = poro[active]
            mesh = extract_submesh(mesh, active)
        else
            poro = max.(poro, minporo)
        end
        return reservoir_domain(mesh, permeability = perm, porosity = poro)
    end

    function check_tiling(tiling)
        if !(tiling isa Tuple{Integer, Integer, Integer}) || any(t -> t < 1, tiling)
            throw(ArgumentError("tiling must be a tuple of three positive integers"))
        end
        return nothing
    end

    function setup_mesh(layers = 1:85; tiling = (1, 1, 1))
        check_tiling(tiling)
        ft = si_unit(:feet)
        nx = 60
        ny = 220
        nz = length(layers)
        cell_dims = (20.0, 10.0, 2.0).*ft
        dims = (nx, ny, nz).*tiling
        mesh = reservoir_mesh(dims, dims.*cell_dims)
        depth = 12000*ft
        offset = [0.0, 0.0, depth]
        for (i, n) in enumerate(mesh.node_points)
            mesh.node_points[i] = n .+ offset
        end
        return mesh
    end

    function spe10_artifact_path()
        return artifact"spe10_model2"
    end

    function setup_system()
        pft3 = si_unit(:pound)/si_unit(:feet)^3
        rhoWS = 64.0*pft3  # water density
        rhoOS = 53.0*pft3  # oil density
        return ImmiscibleSystem(:wo, reference_densities = [rhoWS, rhoOS])
    end

    function setup_wells(domain::Jutul.DataDomain; well_triplets = missing)
        mesh = physical_representation(domain)
        nx, ny, nz = grid_dims_ijk(mesh)
        nx > 0 && nx % 60 == 0 || error("Expected nx to be a positive multiple of 60, got $nx")
        ny > 0 && ny % 220 == 0 || error("Expected ny to be a positive multiple of 220, got $ny")
        if ismissing(well_triplets)
            well_triplets = [
                (30, 110, :I1),
                (1, 1, :P1),
                (60, 1, :P2),
                (60, 220, :P3),
                (1, 220, :P4)
            ]
            tilex, tiley = nx ÷ 60, ny ÷ 220
            ninj = 1
            for j in 0:(tiley-1), i in 0:(tilex-1)
                if i == 0 && j == 0
                    continue
                end
                ninj += 1
                push!(well_triplets, (30 + 60*i, 110 + 220*j, Symbol(:I, ninj)))
            end
            nprod = 4
            for j in 0:tiley, i in 0:tilex
                if i <= 1 && j <= 1
                    continue
                end
                # Keep a single well at each shared corner, using the cell on
                # the lower-index side of internal tile boundaries.
                nprod += 1
                push!(well_triplets, (max(1, 60*i), max(1, 220*j), Symbol(:P, nprod)))
            end
        end
        addw(I, J, name) = setup_vertical_well(domain, I, J, name = name, verbose = false, radius = 5*si_unit(:inch))
        wells = []
        for (I, J, name) in well_triplets
            push!(wells, addw(I, J, name))
        end
        return [w for w in wells]
    end

    function setup_relperm()
        swof = [
            0.200	0.0000	1.0000
            0.250	0.0069	0.8403
            0.300	0.0278	0.6944
            0.350	0.0625	0.5625
            0.400	0.1111	0.4444
            0.450	0.1736	0.3403
            0.500	0.2500	0.2500
            0.550	0.3403	0.1736
            0.600	0.4444	0.1111
            0.650	0.5625	0.0625
            0.700	0.6944	0.0278
            0.750	0.8403	0.0069
            0.800	1.0000	0.0000
        ]
        krw, krow = table_to_relperm(swof, first_label = :w, second_label = :ow)
        return JutulDarcy.ReservoirRelativePermeabilities(w = krw, ow = krow)
    end

    function setup_pvt()
        cP = si_unit(:centi)*si_unit(:poise)
        p_ref = 6000*si_unit(:psi)
        mu_c = 0.0
        mu_ref = 0.3*cP
        b_ref = 1.01
        c_w = 3.10e-6/si_unit(:psi)
        w = JutulDarcy.PVTW(ConstMuBTable(p_ref, b_ref, c_w, mu_ref, mu_c))

        p_o = [300.0, 800.0, 8000.0].*si_unit(:psi)
        b_o = [1.05, 1.02, 1.01]
        mu_o = [2.85, 2.99, 3.0].*cP
        o = JutulDarcy.PVDO(MuBTable(p_o, b_o, mu_o))

        return (w, o)
    end

    function setup_density()
        pvt = setup_pvt()
        return DeckPhaseMassDensities(pvt)
    end

    function setup_viscosity()
        pvt = setup_pvt()
        return DeckPhaseViscosities(pvt)
    end

    function setup_model(;
            layers = 1:85,
            model_arg = NamedTuple(),
            well_triplets = missing,
            domain = missing,
            system = missing,
            kwarg...
        )
        if ismissing(domain)
            domain = setup_reservoir(; layers = layers, kwarg...)
        end
        use_spe10_fluid = ismissing(system)
        if use_spe10_fluid
            system = setup_system()
        end
        wells = setup_wells(domain; well_triplets = well_triplets)
        model = setup_reservoir_model(domain, system; ds_max = 0.1, wells = wells, model_arg...)
        if use_spe10_fluid
            kr = setup_relperm()
            mu = setup_viscosity()
            rho = setup_density()

            rmodel = reservoir_model(model)
            set_secondary_variables!(rmodel,
                RelativePermeabilities = kr,
                PhaseMassDensities = rho,
                PhaseViscosities = mu
            )
            for (k, submodel) in pairs(model.models)
                if JutulDarcy.model_or_domain_is_well(submodel)
                    set_secondary_variables!(submodel,
                        PhaseMassDensities = rho,
                        PhaseViscosities = mu
                    )
                end
            end
        end
        JutulDarcy.set_rock_compressibility!(model, reference_pressure = 6000*si_unit(:psi), compressibility = 1e-6/si_unit(:psi))
        return model
    end

    function setup_state0(model; equil = missing)
        if ismissing(equil)
            p_datum = 6000*si_unit(:psi)
            datum_depth = 12000*si_unit(:feet)
            reservoir = reservoir_domain(model)
            woc = maximum(reservoir[:cell_centroids][3, :]) + 100*si_unit(:feet)
            equil = EquilibriumRegion(model, p_datum, datum_depth, woc = woc)
        elseif equil isa Function
            equil = equil(model)
        end
        return setup_reservoir_state(model, equil)
    end

    """
        setup_well_forces(model; do_scale = true, injection_control = missing)

    Set up controls for all SPE10 wells. Names starting with `I` are injectors;
    other wells are producers. The default injection rate scales with pore volume
    per injector, including vertically repeated layers. A supplied
    `injection_control` is used as-is for every injector.
    """
    function setup_well_forces(model; do_scale = true, injection_control = missing)
        well_names = JutulDarcy.well_symbols(model)
        injectors = filter(name -> startswith(string(name), "I"), well_names)
        if do_scale
            reservoir = reservoir_domain(model)
            pv = pore_volume(reservoir)
            pv_tot = 2.1673386968997316e6
            fraction = sum(pv)/(pv_tot*max(length(injectors), 1))
        else
            fraction = 1.0
        end
        irate = fraction*5000*si_unit(:stb)/si_unit(:day)
        ibhp = 10000*si_unit(:psi)
        pbhp = 4000*si_unit(:psi)

        if ismissing(injection_control)
            sys = reservoir_model(model).system
            ncomp = JutulDarcy.number_of_components(sys)
            mix = zeros(ncomp)
            water = findfirst(isequal(AqueousPhase()), JutulDarcy.get_phases(sys))
            if sys isa MultiPhaseCompositionalSystemLV
                # The immiscible component follows the EOS components.
                component = isnothing(water) ? 1 : ncomp
            else
                component = isnothing(water) ? 1 : water
            end
            mix[component] = 1.0
            phase = isnothing(water) ? 1 : water
            density = JutulDarcy.reference_densities(sys)[phase]
            injection_control = InjectorControl(TotalRateTarget(irate), mix;
                density = density, phases = ((phase, 1.0),))
        end
        pctrl = ProducerControl(BottomHolePressureTarget(pbhp))
        ctrl = Dict{Symbol, Any}()
        limits = Dict{Symbol, Any}()
        for name in well_names
            if name in injectors
                ctrl[name] = injection_control
                limits[name] = (bhp = ibhp, )
            else
                ctrl[name] = pctrl
            end
        end
        return setup_reservoir_forces(model, control = ctrl, limits = limits)
    end
end
