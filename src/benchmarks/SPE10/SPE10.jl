module SPE10
    using LazyArtifacts, Artifacts, JutulDarcy, DelimitedFiles, Jutul

    """
        JutulDarcy.SPE10.setup_case()
        JutulDarcy.SPE10.setup_case(layers = 1:5, nsteps = 200)

    Set up the SPE10, model2, benchmark case. This model is a 1.1 million cell
    nearly incompressible oil-water model with highly heterogenous static
    properties (permeability/porosity).
    """
    function setup_case(;
            num_steps = 100,
            t_total = 2000*si_unit(:day),
            kwarg...
        )
        model = setup_model(; kwarg...)
        state0 = setup_state0(model)
        dt = fill(t_total/num_steps, num_steps)
        forces = setup_well_forces(model)
        forces = setup_well_forces(model)
        return JutulCase(model, dt, forces, state0 = state0)
    end

    function setup_reservoir(;
            layers = 1:85,
            anisotropy = true,
            minporo = 0.01,
            remove_cells = true,
            source = spe10_artifact_path()
        )
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
            return vec(x)
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
        mesh = setup_mesh(layers)
        if remove_cells
            # Remove cells with too low porosity
            active = findall(poro .>= minporo)
            perm = perm[:, active]
            poro = poro[active]
            mesh = extract_submesh(mesh, active)
        else
            poro = max.(poro, minporo)
        end
        return reservoir_domain(mesh, permeability = perm, porosity = poro)
    end

    function setup_mesh(layers = 1:85)
        ft = si_unit(:feet)
        nx = 60
        ny = 220
        nz = length(layers)
        cell_dims = (20.0, 10.0, 2.0).*ft
        dims = (nx, ny, nz)
        mesh = UnstructuredMesh(CartesianMesh(dims, dims.*cell_dims))
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

    function setup_wells(domain::Jutul.DataDomain)
        mesh = physical_representation(domain)
        nx, ny, nz = grid_dims_ijk(mesh)
        nx == 60 || error("Expected nx=60, got $nx")
        ny == 220 || error("Expected ny=220, got $ny")
        addw(I, J, name) = setup_vertical_well(domain, I, J, name = name, verbose = false, radius = 5*si_unit(:inch))
        I1 = addw(30, 110, :I1)
        P1 = addw(1, 1, :P1)
        P2 = addw(60, 1, :P2)
        P3 = addw(60, 220, :P3)
        P4 = addw(1, 220, :P4)
        return [I1, P1, P2, P3, P4]
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

    function setup_model(; layers = 1:85, model_arg = NamedTuple(), kwarg...)
        domain = setup_reservoir(; layers = layers, kwarg...)
        sys = setup_system()
        wells = setup_wells(domain)
        model = setup_reservoir_model(domain, sys; ds_max = 0.1, wells = wells, model_arg...)
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
        JutulDarcy.set_rock_compressibility!(model, reference_pressure = 6000*si_unit(:psi), compressibility = 1e-6/si_unit(:psi))
        return model
    end

    function setup_state0(model)
        p_datum = 6000*si_unit(:psi)
        datum_depth = 12000*si_unit(:feet)
        reservoir = reservoir_domain(model)
        woc = maximum(reservoir[:cell_centroids][3, :]) + 100*si_unit(:feet)
        equil = EquilibriumRegion(model, p_datum, datum_depth, woc = woc)
        return setup_reservoir_state(model, equil)
    end

    function setup_well_forces(model; do_scale = true)
        if do_scale
            reservoir = reservoir_domain(model)
            pv = pore_volume(reservoir)
            pv_tot = 2.1673386968997316e6
            fraction = min(sum(pv)/pv_tot, 1.0)
        else
            fraction = 1.0
        end
        irate = fraction*5000*si_unit(:stb)/si_unit(:day)
        ibhp = 10000*si_unit(:psi)
        pbhp = 4000*si_unit(:psi)

        rhows = JutulDarcy.reference_densities(setup_system())[1]

        ictrl = InjectorControl(TotalRateTarget(irate), [1.0, 0.0], density = rhows)
        pctrl = ProducerControl(BottomHolePressureTarget(pbhp))
        ctrl = Dict(
            :I1 => ictrl,
            :P1 => pctrl,
            :P2 => pctrl,
            :P3 => pctrl,
            :P4 => pctrl,
        )
        limits = Dict{Symbol, Any}(
            :I1 => (bhp = ibhp, )
        )
        return setup_reservoir_forces(model, control = ctrl, limits = limits)
    end
end
