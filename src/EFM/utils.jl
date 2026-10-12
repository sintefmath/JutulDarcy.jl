# Copy data without copying the mesh, which is shared by matrix/fracture domains.
function copy_fracture_data(domain::DataDomain)
    return DataDomain(physical_representation(domain), copy(domain.entities), deepcopy(domain.data))
end

fracture_cell_subset(values::AbstractVector, cells) = values[cells]
fracture_cell_subset(values::AbstractMatrix, cells) = values[:, cells]

function append_cell_property(matrix_values, fracture_values)
    if matrix_values isa AbstractVector && fracture_values isa AbstractVector
        values = vcat(matrix_values, fracture_values)
    else
        mv = matrix_values
        fv = fracture_values
        if mv isa AbstractVector
            mv = reshape(mv, 1, :)
        end
        if fv isa AbstractVector
            fv = reshape(fv, 1, :)
        end
        if size(mv, 1) != size(fv, 1)
            # Isotropic fracture properties in the matrix's symmetric tensor format.
            if size(fv, 1) == 1 && size(mv, 1) in (3, 6)
                fv = repeat(fv, size(mv, 1), 1)
                if size(mv, 1) == 6
                    fv[4:6, :] .= 0
                end
            else
                throw(ArgumentError("matrix and fracture cell property dimensions do not match"))
            end
        end
        values = hcat(mv, fv)
    end
    return values
end

# Property overrides have identical meaning for EFM patches and cut-cell DFM.
function fracture_properties(objects, ids, ::Type{T}; aperture = nothing,
        hydraulic_aperture = nothing, permeability = nothing, porosity = nothing) where T
    defaults = fracture_property_defaults(objects, T)
    a = fracture_property(aperture, defaults.aperture[ids], ids, "aperture")
    k = fracture_property(permeability, defaults.permeability[ids], ids, "permeability")
    phi = fracture_property(porosity, defaults.porosity[ids], ids, "porosity")
    ah = a
    if !isnothing(hydraulic_aperture)
        ah = fracture_property(hydraulic_aperture, a, ids, "hydraulic_aperture")
        all(>=(0), ah) || throw(ArgumentError("hydraulic_aperture must be non-negative"))
    end
    if isnothing(permeability) && (!isnothing(hydraulic_aperture) || !isnothing(aperture))
        k = ah.^2/12
    end
    all(>(0), a) || throw(ArgumentError("aperture must be positive"))
    all(>=(0), k) || throw(ArgumentError("permeability must be non-negative"))
    all(x -> 0 < x <= 1, phi) || throw(ArgumentError("porosity must lie in (0, 1]"))
    return (aperture = a, permeability = k, porosity = phi)
end

function embedded_cell_faces(neighbors::AbstractMatrix, nc)
    faces, pos = Jutul.get_facepos(neighbors, nc)
    return Jutul.IndirectionMap(faces, pos)
end

function fracture_patch_spatial_index(patches::Vector{FracturePatch{T}}, frames) where T
    low, high = Point3{T}[], Point3{T}[]
    for patch in patches
        frame = frames[patch.fracture]
        points = [unproject(frame, p) for part in patch.parts for p in part]
        push!(low, reduce((a, b) -> min.(a, b), points))
        push!(high, reduce((a, b) -> max.(a, b), points))
    end
    return CellSpatialIndex(low, high)
end
