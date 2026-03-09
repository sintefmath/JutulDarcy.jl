function Jutul.triangulate_mesh(m::DataDomain{T}; kwarg...) where {T<:SimpleWell}
    return triangulate_well_domain(m; is_simple = true, kwarg...)
end

function Jutul.triangulate_mesh(m::DataDomain{T}; kwarg...) where {T<:MultiSegmentWell}
    return triangulate_well_domain(m; is_simple = false, kwarg...)
end

function triangulate_well_domain(w::DataDomain;
        rfactor = 10.0,
        outer = false,
        is_simple = false,
        flatten = true,
        flip = false, type = Cells()
    )
    is_perf = type == JutulDarcy.Perforations()
    if is_simple
        global_pts = w[:perforation_centroids]
        branches = [1:size(global_pts, 2)]
    else
        global_pts = w[:cell_centroids]
        N = w.representation.neighborship
        branches = find_well_branches(N, overlap = true)
    end
    # if size(pts, 2) == 1
    #     pts = repeat(pts, 1, 2)
    #     pt_remap = [1, 1]
    # else
    #     pt_remap = 1:size(pts, 2)
    # end

    offset = 0
    allpts = []
    alltri = []
    allcells = Int[]
    for (branch_no, branch) in enumerate(branches)
        pts = global_pts[:, branch]
        if is_simple && !is_perf
            pt_remap = fill(1, size(pts, 2))
        else
            pt_remap = collect(1:size(pts, 2))
        end
        if size(pts, 1) == 3 && branch_no == 1
            pushfirst!(pt_remap, 1)
            pts = hcat(pts[:, 1], pts)
            pts[3, 1] -= 10.0
        end

        radius = w[:radius] .* rfactor
        r = sum(radius) / length(radius)

        ctri = triangulate_cylindrical_mesh(pts, radius = r, cap = :round)
        cell_index = pt_remap[ctri.vert_srcidx]
        if is_simple && !is_perf
            cell_index .= 1
        end
        push!(allpts, ctri.xyz')
        push!(alltri, ctri.tri' .+ offset)
        append!(allcells, branch[cell_index])
        offset += size(ctri.xyz, 2)

    end
    allpts = vcat(allpts...)
    alltri = vcat(alltri...)

    mapper = (
        Cells = (cell_data) -> cell_data[allcells],
        indices = (Cells = allcells, )
    )
    return (points = allpts, triangulation = alltri, mapper = mapper)
end

function find_well_branches(N; overlap = false)
    N = reinterpret(Tuple{Int, Int}, N)
    Nmax = maximum(t -> maximum(t), N)
    counts = zeros(Int, Nmax)
    for (l, r) in N
        counts[l] += 1
        counts[r] += 1
    end
    terminus = findall(x -> x == 1, counts)

    function dfs_color(N, start, visited)
        branch = Int[]
        return dfs_color!(branch, N, start, visited)
    end

    function dfs_color!(branch, N, start, visited)
        visited[start] = true
        push!(branch, start)
        conn_idx = findfirst(t -> start in t, N)
        l, r = N[conn_idx]
        if l == start && !visited[r]
            dfs_color!(branch, N, r, visited)
        elseif r == start && !visited[l]
            dfs_color!(branch, N, l, visited)
        end
        return branch
    end
    visited = falses(Nmax)
    branches = Vector{Vector{Int}}()

    while length(terminus)> 0
        start = popfirst!(terminus)
        new_branch = dfs_color(N, start, visited)
        if length(branches) != 0
            reverse!(new_branch)
        end
        push!(branches, new_branch)
    end
    if overlap
        sort_tup(a) = ifelse(a[1] > a[2], (a[2], a[1]), a)
        for (bno, branch) in enumerate(branches)
            for bother in branches
                if bother == branch
                    continue
                end
                start = branch[1]
                self = sort_tup((start, bother[end]))
                if !isnothing(self)
                    if self[1] == start
                        other = self[2]
                    else
                        other = self[1]
                    end
                    pushfirst!(branch, other)
                end
            end
        end
    end
    @info branches N
    return branches
end
