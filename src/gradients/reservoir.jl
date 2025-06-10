
function setup_reservoir_parameter_optimization_generic(case::JutulCase;
        use_trans = false,
        use_pore_volume = false,
        strict = false,
        kwarg...
    )
    is_multimodel = case.model isa MultiModel
    if is_multimodel
        rparameters = case.parameters[:Reservoir]
    else
        rparameters = case.parameters
    end
    rmodel = reservoir_model(case.model)
    rdomain = reservoir_domain(rmodel)
    # TODO: Rename...
    # data_domain (non-geometry)
    # wells - WI

    skip_list = [
        :satnum,
        :eqlnum,
        :pvtnum,
        :cell_centroids,
        :volumes,
        :neighbors,
        :areas,
        :normals,
        :face_centroids,
        :half_face_cells,
        :half_face_faces,
        :boundary_areas,
        :boundary_centroids,
        :boundary_normals,
        :boundary_neighbors,
    ]

    DT = OrderedDict{Symbol, Any}
    opt_dict = DT()
    if use_trans
        opt_dict[:transmissibilities] = copy(rparameters[:Transmissibilities])
        push!(skip_list, :permeability)
    end
    if use_pore_volume
        opt_dict[:pore_volume] = pore_volume(rdomain)
        push!(skip_list, :porosity)
        push!(skip_list, :net_to_gross)
    end
    for k in setdiff(keys(rdomain), skip_list)
        v = rdomain[k]
        if eltype(v)<:AbstractFloat
            opt_dict[k] = v
        end
    end
    if is_multimodel
        for (k, submodel) in pairs(case.model.models)
            if model_or_domain_is_well(submodel)
                subdict = DT()
                subprm = case.parameters[k]
                subdict[:WellIndices] = subprm[:WellIndices]
                if haskey(subprm, :WellIndicesThermal)
                    subdict[:WellIndicesThermal] = subprm[:WellIndicesThermal]
                end
                opt_dict[k] = subdict
            end
        end
    end
    dp = DictParameters(opt_dict, strict = strict)
    F = missing
    return (dp, F)
end

