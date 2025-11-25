module SPE10
    using LazyArtifacts, Artifacts, JutulDarcy, DelimitedFiles, Jutul
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
        md = si_unit(:darcy)

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
            offset = nvals_data ÷ 6
            vals = perm[:, ((i-1)*offset + 1):(i*offset)]
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
end
