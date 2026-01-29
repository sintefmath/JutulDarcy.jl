function setup_well(g_matrix::JutulMesh, K_matrix, matrix_cells::AbstractVector,
    g_fractures::JutulMesh, K_fractures, matrix_fracture_faces;
    simple_well = false,
    kwarg...
    )

    simple_well == false || throw(ArgumentError("Only MultiSegmentWell setup is supported for DFM wells."))
    tmp = setup_well(g_matrix, K_matrix, matrix_cells; simple_well=simple_well, kwarg...)

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
    W = MultiSegmentWell(
            neighborship,
            tmp.representation.perforations.reservoir, 
            tmp.representation.perforations.self,
            perforation_cells_fractures,
            perforation_cells_fractures_self;
            end_nodes = tmp.representation.end_nodes,
            type = tmp.representation.type,
            name = tmp.representation.name,
            segment_models = segment_models,
            surface_conditions = tmp.representation.surface,
    )

    Wdomain = DataDomain(W)
    c = Cells()
    p = Perforations()
    f = Faces()

    

    return W, tmp

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
        perforation_cells_matrix, perforation_cells_matrix_self,
        perforation_cells_fractures, perforation_cells_fractures_self;
        end_nodes = missing,
        type = :ms,
        name = :Well,
        segment_models = nothing,
        surface_conditions = default_surface_cond(),
    )
    size(neighbors, 1) == 2 || throw(ArgumentError("Connectivity matrix for multisegment well must have two rows"))
    num_nodes = maximum(neighbors, init = 1)
    num_segments = size(neighbors, 2)
    num_perf = length(perforation_cells_matrix_self) + length(perforation_cells_fractures_self)
    maximum(vcat(perforation_cells_matrix_self, perforation_cells_fractures_self)) <= num_nodes || throw(ArgumentError("Perforation cells (self) must be less than or equal to number of nodes $(num_nodes), was $(maximum(vcat(perforation_cells_matrix_self, perforation_cells_fractures_self)))"))
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
    perf = (self = perforation_cells_matrix_self, reservoir = perforation_cells_matrix,
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

