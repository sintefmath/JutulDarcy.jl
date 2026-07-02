
function JutulDarcy.plot_faults!(ax, mesh::UnstructuredMesh; domain = missing, textcolor = :black, fontsize = 12, toggle = missing, kwarg...)
    faults = get_mesh_entity_tag(mesh, Faces(), :faults, throw = false)
    if !ismissing(faults)
        n = length(keys(faults))
        i = 1
        for (k, v) in faults
            if length(v) == 0
                continue
            end
            flt = plot_mesh!(ax, mesh; faces = v, color = i, colorrange = (1, max(n, 2)), kwarg...)
            if !ismissing(domain) && fontsize > 0
                center = domain[:face_centroids][:, v[1]]
                for i in 2:length(v)
                    center += domain[:face_centroids][:, v[i]]
                end
                center /= length(v)
                txt = text!(ax, center...; text = "$k", color = textcolor, fontsize = fontsize)
            else
                txt = missing
            end
            if !ismissing(toggle)
                on(toggle.checked) do checked
                    flt.visible[] = checked
                    if !ismissing(txt)
                        txt.visible[] = checked
                    end
                end
            end
            i += 1
        end
    end
    return ax
end
