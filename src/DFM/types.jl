struct FracturePerforations <: JutulEntity end
struct FractureMatrixConnection <: JutulEntity end

abstract type AbstractFractureMatrixConductivity <: ScalarVariable end

struct FractureMatrixTransmissibility <: AbstractFractureMatrixConductivity end
Jutul.variable_scale(::FractureMatrixTransmissibility) = 1e-10
Jutul.minimum_value(::FractureMatrixTransmissibility) = 0.0
Jutul.associated_entity(::FractureMatrixTransmissibility) = FractureMatrixConnection()

function Jutul.default_parameter_values(data_domain, model, param::FractureMatrixTransmissibility, symb)

    fmc = FractureMatrixConnection()
    if haskey(data_domain, :fracture_matrix_transmissibility, fmc)
        return data_domain[:fracture_matrix_transmissibility, fmc]
    end

    # Get matrix-fracture connection target/source cells
    fracture_cells = data_domain[:connection_cells, fmc]
    # Sanity check
    fmesh = physical_representation(data_domain)
    if hasproperty(fmesh, :intersection_cells)
       .!any(c ∈ fmesh.intersection_cells for c in fracture_cells) || error(
        "Fracture cells must not include intersection cells for matrix-fracture\
        cross term setup.")
    end
    # Compute transmissibilities
    T = matrix_fracture_connection_conductivity(data_domain, :permeability)
    return T

  end

struct FractureMatrixThermalConductivity <: AbstractFractureMatrixConductivity end
Jutul.variable_scale(::FractureMatrixThermalConductivity) = 1.0
Jutul.minimum_value(::FractureMatrixThermalConductivity) = 0.0
Jutul.associated_entity(::FractureMatrixThermalConductivity) = FractureMatrixConnection()

function Jutul.default_parameter_values(data_domain, model, param::FractureMatrixThermalConductivity, symb)
    
    fmc = FractureMatrixConnection()
    if haskey(data_domain, :fracture_matrix_thermal_conductivity, fmc)
        return data_domain[:fracture_matrix_thermal_conductivity, fmc]
    end

    # Get matrix-fracture connection target/source cells
    fracture_cells = data_domain[:connection_cells, fmc]
    # Sanity check
    fmesh = physical_representation(data_domain)
    if hasproperty(fmesh, :intersection_cells)
       .!any(c ∈ fmesh.intersection_cells for c in fracture_cells) || error(
        "Fracture cells must not include intersection cells for matrix-fracture\
        cross term setup.")
    end
    # Compute transmissibilities
    T = matrix_fracture_connection_conductivity(data_domain, :thermal_conductivity)
    return T

end

struct FractureMatrixGravityDifference <: ScalarVariable end
Jutul.variable_scale(::FractureMatrixGravityDifference) = 1.0
Jutul.minimum_value(::FractureMatrixGravityDifference) = -Inf
Jutul.associated_entity(::FractureMatrixGravityDifference) = FractureMatrixConnection()

function Jutul.default_parameter_values(data_domain, model, param::FractureMatrixGravityDifference, symb)
    
    fmc = FractureMatrixConnection()
    if haskey(data_domain, :fracture_matrix_gravity_difference, fmc)
        return data_domain[:fracture_matrix_gravity_difference, fmc]
    end

    # Get matrix-fracture connection target/source cells
    fracture_cells = data_domain[:connection_cells, fmc]
    # Sanity check
    fmesh = physical_representation(data_domain)
    if hasproperty(fmesh, :intersection_cells)
       .!any(c ∈ fmesh.intersection_cells for c in fracture_cells) || error(
        "Fracture cells must not include intersection cells for matrix-fracture\
        cross term setup.")
    end
    # Compute gravity potential difference
    xm = data_domain[:matrix_cell_centroids, fmc]
    xf = data_domain[:cell_centroids, Jutul.Cells()][:, fracture_cells]
    g = gravity_constant
    gdz = g * (xm[3, :] .- xf[3, :])

    return gdz

end

struct FractureWellIndices <: ScalarVariable end

Jutul.minimum_value(::FractureWellIndices) = 0.0
Jutul.variable_scale(::FractureWellIndices) = 1e-10
Jutul.associated_entity(::FractureWellIndices) = FracturePerforations()

function Jutul.default_parameter_values(data_domain, model, param::FractureWellIndices, symb)
    if count_entities(data_domain, FracturePerforations()) == 0
        return Float64[]
    end
    WI = copy(data_domain[:well_index_frac, FracturePerforations()])
    dims = data_domain[:cell_dims_frac, FracturePerforations()]
    perm = data_domain[:permeability_frac, FracturePerforations()]
    net_to_gross = data_domain[:net_to_gross_frac, FracturePerforations()]
    direction = data_domain[:perforation_direction_frac, FracturePerforations()]
    skin = data_domain[:skin_frac, FracturePerforations()]
    Kh = data_domain[:Kh_frac, FracturePerforations()]
    radius = data_domain[:perforation_radius_frac, FracturePerforations()]
    drainage_radius = data_domain[:drainage_radius_frac, FracturePerforations()]
    gdim = size(data_domain[:cell_centroids, Cells()], 1)

    T = Base.promote_type(eltype(WI), eltype(skin), eltype(radius),
    eltype(perm), eltype(Kh), eltype(net_to_gross), eltype(drainage_radius))
    if T != eltype(WI)
        WI = convert(Vector{T}, WI)
    end

    for (i, val) in enumerate(WI)
        defaulted = !isfinite(val)
        if defaulted
            Δ = dims[i]
            if perm isa AbstractVector
                K = perm[i]
            else
                K = perm[:, i]
            end
            K = Jutul.expand_perm(K, gdim)
            r = radius[i]
            dir = direction[i]
            WI[i] = compute_peaceman_index(Δ, K, r, dir;
                skin = skin[i],
                Kh = Kh[i],
                net_to_gross = net_to_gross[i],
                drainage_radius = drainage_radius[i]
            )
        end
    end
    return WI
end

struct FractureWellIndicesThermal <: ScalarVariable end

Jutul.minimum_value(::FractureWellIndicesThermal) = 0.0
Jutul.variable_scale(::FractureWellIndicesThermal) = 1.0
Jutul.associated_entity(::FractureWellIndicesThermal) = FracturePerforations()
function Jutul.default_parameter_values(data_domain, model, param::FractureWellIndicesThermal, symb)
    WIt = copy(data_domain[:thermal_well_index_frac, FracturePerforations()])
    dims = data_domain[:cell_dims_frac, FracturePerforations()]
    thermal_conductivity = data_domain[:thermal_conductivity_frac, FracturePerforations()]
    direction = data_domain[:perforation_direction_frac, FracturePerforations()]
    radius = data_domain[:perforation_radius_frac, FracturePerforations()]
    drainage_radius = data_domain[:drainage_radius_frac, FracturePerforations()]
    gdim = size(data_domain[:cell_centroids, Cells()], 1)

    # These are defined per cell, map to perforations
    well = physical_representation(data_domain)
    ic = well.perforations.self_fracture

    λ_casing = data_domain[:casing_thermal_conductivity, Cells()][ic]
    λ_grout = data_domain[:grouting_thermal_conductivity, Cells()][ic]
    casing_thickness = data_domain[:casing_thickness, Cells()][ic]
    grouting_thickness = data_domain[:grouting_thickness, Cells()][ic]

    T = Base.promote_type(
        eltype(WIt), eltype(thermal_conductivity), eltype(radius),
        eltype(λ_casing), eltype(λ_grout), eltype(casing_thickness),
        eltype(grouting_thickness), eltype(drainage_radius))
    if T != eltype(WIt)
        WIt = convert(Vector{T}, WIt)
    end

    for (i, val) in enumerate(WIt)
        defaulted = !isfinite(val)
        if defaulted
            Δ = dims[i]
            if thermal_conductivity isa AbstractVector
                Λ_i = thermal_conductivity[i]
            else
                Λ_i = thermal_conductivity[:, i]
            end
            Λ_i = Jutul.expand_perm(Λ_i, gdim)
            WIt[i] = compute_well_thermal_index(Δ, Λ_i, radius[i], direction[i];
                casing_thickness = casing_thickness[i],
                grouting_thickness = grouting_thickness[i],
                casing_thermal_conductivity = λ_casing[i],
                grouting_thermal_conductivity = λ_grout[i],
                drainage_radius = drainage_radius[i],
            )
        end
    end
    return WIt
end