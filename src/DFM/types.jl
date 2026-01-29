struct FracturePerforations <: JutulEntity end

struct FractureWellIndices <: ScalarVariable end

Jutul.minimum_value(::FractureWellIndices) = 0.0
Jutul.variable_scale(::FractureWellIndices) = 1e-10

Jutul.associated_entity(::FractureWellIndices) = FracturePerforations()

function Jutul.default_parameter_values(data_domain, model, param::FractureWellIndices, symb)
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
    for (i, val) in enumerate(WI)
        defaulted = !isfinite(val)
        if defaulted
            error("Fracture well indices cannot be defaulted")
        end
    end
    return WI
end