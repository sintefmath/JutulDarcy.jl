function setup_well(g_matrix::JutulMesh, K_matrix, matrix_cells::AbstractVector,
    g_fractures::JutulMesh, K_fractures, matrix_fracture_faces;
    simple_well = false,
    frac_args = NamedTuple(),
    kwarg...
    )

    !simple_well || throw(ArgumentError("Only MultiSegmentWell setup is supported for DFM wells."))
    tmp = setup_well(g_matrix, K_matrix, matrix_cells; simple_well=simple_well, kwarg...)

    # return tmp

    neighborship, perforation_cells_fractures, perforation_cells_fractures_self, old_ix = 
    add_fracture_cells(g_matrix, tmp.representation, matrix_fracture_faces)

    # return neighborship, perforation_cells_fractures, perforation_cells_fractures_self, old_ix, tmp

    num_segments = size(neighborship, 2)
    segment_models = [SegmentWellBoreFrictionHB() for _ in 1:num_segments]
    for seg = 1:num_segments
        if old_ix[seg] != 0
            segment_models[seg] = tmp.representation.segment_models[old_ix[seg]]
        end
    end

    # perforation_cells = vcat(
    #     tmp.representation.perforations.reservoir,
    #     perforation_cells_fractures
    # )
    # perforation_cells_self = vcat(
    #     tmp.representation.perforations.self,
    #     perforation_cells_fractures_self
    # )
    
    # nrp = length(tmp.representation.perforations.reservoir)
    # fracture_perforations = collect(nrp+1:length(perforation_cells))
    perf = (
        self = tmp.representation.perforations.self,
        reservoir = tmp.representation.perforations.reservoir,
        self_fracture = perforation_cells_fractures_self,
        fracture = perforation_cells_fractures,
    )

    tmp_frac = setup_well(g_fractures, K_fractures, perforation_cells_fractures;
    simple_well=simple_well, frac_args...)

    W = MultiSegmentWell(
        tmp.representation.type,
        maximum(neighborship, init = 1),
        length(segment_models),
        length(tmp.representation.perforations.self),
        perf,
        neighborship,
        tmp.representation.end_nodes,
        tmp.representation.surface,
        tmp.representation.name,
        segment_models,
    )

    # TODO: Make FracturePerforations() entity type and move all
    # fracture-related perforation logic there + make FractureFromWellFlowCT

    Wdomain = DataDomain(W)
    
    for (k, v) in tmp.data
        val, entity = v
        frac_val = tmp_frac[k, entity]
        if entity == Cells()
            # Cell data
            if k ∈ [:cell_centroids, :perforation_centroids]
                new_val = hcat(val, frac_val)
            else
                new_val = vcat(val, frac_val)
            end
        elseif entity == Faces()
            new_val = similar(val, count_entities(W, Faces()))
            for seg in 1:num_segments
                if old_ix[seg] != 0
                    new_val[seg] = val[old_ix[seg]]
                else
                    if old_ix[seg-1] > 0
                        new_val[seg] = val[old_ix[seg-1]]
                    elseif old_ix[seg+1] > 0
                        new_val[seg] = val[old_ix[seg+1]]
                    end
                end
            end
        else
            new_val = val
        end
        Wdomain[k, entity] = new_val
    end

    return Wdomain

end

function add_fracture_cells(g_matrix, representation, fracture_faces)
    N = representation.neighborship
    nc = maximum(N)
    w2m = fill(0, nc)
    perf = representation.perforations
    w2m[perf.self] .= perf.reservoir

    N_new = Matrix{Int}(undef, 2, 0)
    N_matrix = get_neighborship(g_matrix)
    perforation_cells_fractures = Int[]
    perforation_cells_fractures_self = Int[]

    wc = maximum(N)
    added = false
    old_ix = Int[]
    for (seg, mcells) in enumerate(eachcol(w2m[N]))
        ix = vec(all(mcells .== N_matrix, dims = 1) .|| all(reverse(mcells) .== N_matrix, dims = 1))
        if any(ix)
            @assert sum(ix) == 1 "Failed to find unique matching face for well segment perforation"
            face = findfirst(ix)
            if face ∈ fracture_faces
                wc += 1
                fcell = findfirst(fracture_faces .== face)
                push!(perforation_cells_fractures, fcell)
                push!(perforation_cells_fractures_self, wc)
                m2f = [N[1, seg], wc]
                f2m = [wc, N[2, seg]]
                N_new = hcat(N_new, m2f, f2m)
                push!(old_ix, 0, 0)
                added = true
            end
        end
        if !added
            N_new = hcat(N_new, N[:, seg])
            push!(old_ix, seg)
        end
        added = false
    end

    return N_new, perforation_cells_fractures, perforation_cells_fractures_self, old_ix
end

function MultiSegmentWell(neighbors::AbstractMatrix,
        perforation_cells_reservoir, perforation_cells_self, fracture_perforations;
        end_nodes = missing,
        type = :ms,
        name = :Well,
        segment_models = nothing,
        surface_conditions = default_surface_cond(),
    )
    size(neighbors, 1) == 2 || throw(ArgumentError("Connectivity matrix for multisegment well must have two rows"))
    num_nodes = maximum(neighbors, init = 1)
    num_segments = size(neighbors, 2)
    num_perf = length(perforation_cells_self) + length(fracture_perforations)
    maximum(vcat(perforation_cells_self, fracture_perforations)) <= num_nodes || throw(ArgumentError("Perforation cells (self) must be less than or equal to number of nodes $(num_nodes), was $(maximum(vcat(perforation_cells_self, fracture_perforations)))"))
    if ismissing(end_nodes)
        from_nodes = unique(neighbors[1, :])
        to_nodes = unique(neighbors[2,:])
        end_nodes = setdiff(to_nodes, from_nodes)
        if isempty(end_nodes)
            @assert num_nodes == 1 "Failed to determine end_nodes from connectivity matrix, please provide explicitly."
            end_nodes = [1]
        end
    end
    if isnothing(segment_models)
        segment_models = [SegmentWellBoreFrictionHB() for _ in 1:num_segments]
    else
        segment_models::AbstractVector
        length(segment_models) == num_segments || throw(ArgumentError("If segment models are provided, there must be one per segment. Was $(length(segment_models)), expected $num_segments"))
    end
    perf = (self = perforation_cells_self, reservoir = perforation_cells_reservoir,
            self_fracture = perforation_cells_fractures_self, fracture = perforation_cells_fractures)
    return MultiSegmentWell(
        type,
        num_nodes,
        num_segments,
        num_perf,
        perf,
        neighbors,
        end_nodes,
        surface_conditions,
        name,
        segment_models,
    )
end
