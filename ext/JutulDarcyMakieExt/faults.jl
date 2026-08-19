
function JutulDarcy.plot_faults!(ax, mesh::JutulMesh; domain = missing, textcolor = :black, fontsize = 12, toggle = missing, kwarg...)
    faults = JutulDarcy.get_faults(mesh)
    is_coarse = mesh isa CoarseMesh
    if !ismissing(faults)
        n = length(keys(faults))
        i = 1
        for (k, v) in faults
            if length(v) == 0
                continue
            end
            if is_coarse
                coarse_faces = v
                fault_faces = Int[]
                for cf in coarse_faces
                    append!(fault_faces, mesh.coarse_faces_to_fine[cf])
                end
                plot_mesh = mesh.parent
            else
                plot_mesh = mesh
                fault_faces = v
            end
            flt = plot_mesh!(ax, plot_mesh;
                faces = fault_faces,
                color = i,
                colorrange = (1, max(n, 2)),
                kwarg...
            )
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
