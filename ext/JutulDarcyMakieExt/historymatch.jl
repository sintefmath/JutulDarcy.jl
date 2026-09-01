function JutulDarcy.plot_mismatch(obj, res::ReservoirSimResult)
    _, w = JutulDarcy.HistoryMatching.evaluate_match(obj, res, log = true)
    return JutulDarcy.plot_mismatch(w)
end

function JutulDarcy.plot_mismatch(w)
    (; impact, targets, wells) = JutulDarcy.HistoryMatching.compute_well_target_contribution_matrix(w)
    nw, ntargets = size(impact)

    fig = Figure(size = (1400, 900))
    # Barplot over well contributions, banked by control type
    ax1 = Axis(fig[2, 1], title = "Contribution by well (total)")
    ax2 = Axis(fig[2, 2], title = "Contribution by target (total)")
    ax3 = Axis(fig[3, 1:2], title = "Contribution by target (per step)")
    ax4 = Axis(fig[4, 1:2], title = "Contribution by well (per step)")

    colors = Makie.wong_colors()

    well_idx = repeat((1:nw), 1, ntargets)
    target_idx = repeat((1:ntargets)', nw, 1)
    tvec = vec(target_idx)
    wvec = vec(well_idx)
    barplot!(ax1, wvec, vec(impact), stack = wvec, color = colors[tvec])

    labels = map(String, targets)
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    Legend(fig[1, 1:2], elements, labels, orientation  = :horizontal, tellheight = true)
    ax1.xticks[] = (1:nw, map(string, wells))
    ax1.xticklabelrotation[] = -pi/4

    # Barplot over control types
    target_weight = sum(impact, dims = 1)
    barplot!(ax2, 1:ntargets, vec(target_weight), color = colors[1:ntargets])
    ax2.xticks[] = (1:ntargets, map(string, targets))

    # Magnitude of control by step
    maxstep = 0
    for well in wells
        for wl in w[well]
            maxstep = max(maxstep, wl.stop)
        end
    end
    target_step_vals = zeros(ntargets, maxstep)
    well_step_vals = zeros(nw, maxstep)
    for step in 1:maxstep
        (; impact, targets, wells) = JutulDarcy.HistoryMatching.compute_well_target_contribution_matrix(w, step = step)
        target_step_vals[:, step] = sum(impact, dims = 1)'
        well_step_vals[:, step] = sum(impact, dims = 2)
    end
    for i in 1:ntargets
        lines!(ax3, 1:maxstep, target_step_vals[i, :], color = colors[i], label = string(targets[i]))
    end
    totals = vec(sum(target_step_vals, dims = 1))
    lines!(ax3, 1:maxstep, totals, color = :black, label = "Total")
    axislegend(ax3)
    # Magnitude of well by step
    plts = []
    for i in 1:nw
        plt = lines!(ax4, 1:maxstep, well_step_vals[i, :], label = string(wells[i]))
        push!(plts, plt)
    end
    hidexdecorations!(ax3)
    Legend(fig[5, 1:2], plts, map(string, wells), orientation  = :horizontal, tellheight = true, nbanks = 3)

    fig
end
