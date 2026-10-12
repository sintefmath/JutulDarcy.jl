"""
    add_thermal_dispersion!(domain::DataDomain, D)

Add thermal dispersivity to a reservoir domain. `D` is either a 2×1 matrix
`[α_L; α_T]` (longitudinal and transverse dispersivity, used for all cells) or
a 1×nc matrix with a single dispersivity per cell (used for both components).

The dispersive heat flux is `-ρ_w c_w D ∇T` with
`D = α_T |v| I + (α_L - α_T) v v'/|v|`, where `v` is the total Darcy velocity.
`D` is evaluated explicitly at the start of each time step from the velocity
reconstructed from the face fluxes, and enters the energy equation as an
additional face transmissibility.
"""
function add_thermal_dispersion!(domain::DataDomain, D::AbstractArray{<:Real, 2})
    D = if size(D, 1) == 2 && size(D, 2) == 1
        repeat(D, 1, number_of_cells(domain))
    elseif size(D, 1) == 1 && size(D, 2) == number_of_cells(domain)
        repeat(D, 2, 1)
    else
        error("Unsupported size for thermal dispersivity array")
    end
    domain[:thermal_dispersivity] = D
    return domain
end

# Dispersion is computed in three stages before each time step, in the order
# the parameters are added to the model:
#   1. TotalVolumeFlux: total Darcy volume flux over each face
#   2. TotalDarcyVelocity: cell velocity reconstructed from the face fluxes
#   3. ThermalDispersionTransmissibility: half-face transmissibilities of D
# All topology is taken from the model's own discretization and all geometry
# from the face parameter ThermalDispersionGeometry, so that the computations
# only use data that is local to the model.

"""
    ThermalDispersionGeometry(dim)

Static face geometry for thermal dispersion. For each face, rows `1:dim` and
`dim+1:2dim` hold the vectors from the left and right cell centroids to the face
centroid, and rows `2dim+1:3dim` the area-weighted normal (oriented from left to
right).
"""
struct ThermalDispersionGeometry{D} <: VectorVariables end
ThermalDispersionGeometry(dim::Int) = ThermalDispersionGeometry{dim}()
Jutul.values_per_entity(model, ::ThermalDispersionGeometry{D}) where D = 3*D
Jutul.associated_entity(::ThermalDispersionGeometry) = Faces()

function Jutul.default_parameter_values(data_domain, model, ::ThermalDispersionGeometry{D}, symb) where D
    N = data_domain[:neighbors]
    xc = data_domain[:cell_centroids]
    xf = data_domain[:face_centroids]
    normals = data_domain[:normals]
    areas = data_domain[:areas]
    nf = size(N, 2)
    G = zeros(3*D, nf)
    for face in 1:nf
        l, r = N[1, face], N[2, face]
        for i in 1:D
            G[i, face] = xf[i, face] - xc[i, l]
            G[D + i, face] = xf[i, face] - xc[i, r]
            G[2D + i, face] = areas[face]*normals[i, face]
        end
    end
    return G
end

@inline function dispersion_face_geometry(G, face, ::Val{D}) where D
    C_l = SVector{D, Float64}(ntuple(i -> @inbounds(G[i, face]), Val(D)))
    C_r = SVector{D, Float64}(ntuple(i -> @inbounds(G[D + i, face]), Val(D)))
    NA = SVector{D, Float64}(ntuple(i -> @inbounds(G[2D + i, face]), Val(D)))
    return (C_l, C_r, NA)
end

dispersion_neighborship(g::MinimalTPFATopology) = g.neighborship
dispersion_neighborship(g) = get_neighborship(g)

"""
Cell-to-face connections `(conn_pos, conn_data)` of the model, where
`conn_data[conn_pos[c]:conn_pos[c+1]-1]` are the faces of cell `c` with entries
`(self, other, face, face_sign)`.
"""
function dispersion_cell_connections(model)
    discr = model.domain.discretizations
    disc = haskey(discr, :heat_flow) ? discr.heat_flow : discr.mass_flow
    if !(disc isa TwoPointPotentialFlowHardCoded)
        N = dispersion_neighborship(physical_representation(model.domain))
        disc = TwoPointPotentialFlowHardCoded(N, number_of_cells(model.domain))
    end
    return (disc.conn_pos, disc.conn_data)
end

"""
    TotalVolumeFlux()

Total Darcy volume flux over each face (sum over phases), evaluated at the start
of each time step.
"""
struct TotalVolumeFlux <: ScalarVariable end
Jutul.associated_entity(::TotalVolumeFlux) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::TotalVolumeFlux, symb)
    return zeros(number_of_faces(data_domain))
end

function Jutul.update_parameter_before_step!(q_t, ::TotalVolumeFlux, storage, model, dt, forces)
    N = dispersion_neighborship(physical_representation(model.domain))
    update_total_volume_flux!(q_t, N, storage.state, model)
    return q_t
end

function update_total_volume_flux!(q_t, N, state, model)
    F = face -> @inbounds q_t[face] = face_total_volume_flux(face, N, state, model)
    Jutul.threaded_loop_minbatch(F, length(q_t), model.context)
    return q_t
end

@inline function face_total_volume_flux(face, N, state, model)
    @inbounds l, r = N[1, face], N[2, face]
    # TODO: This assumes the default discretizations.
    tpfa = TPFA(l, r, 1)
    upw = SPU(l, r)
    q_face = darcy_phase_volume_fluxes(face, state, model, Jutul.DefaultFlux(), tpfa, upw)
    q_sum = 0.0
    for q in q_face
        q_sum += value(q)
    end
    return q_sum
end

"""
    TotalDarcyVelocity(dim)

Cell-wise total Darcy velocity reconstructed from the face fluxes as
`v_c = Σ_f q_f (x_f - x_c) / V_c`, evaluated at the start of each time step.
Fluxes over boundary faces and well connections are not included.
"""
struct TotalDarcyVelocity{D} <: VectorVariables end
TotalDarcyVelocity(dim::Int) = TotalDarcyVelocity{dim}()
Jutul.values_per_entity(model, ::TotalDarcyVelocity{D}) where D = D
Jutul.associated_entity(::TotalDarcyVelocity) = Cells()

function Jutul.default_parameter_values(data_domain, model, param::TotalDarcyVelocity{D}, symb) where D
    return zeros(D, number_of_cells(data_domain))
end

function Jutul.update_parameter_before_step!(v, ::TotalDarcyVelocity{D}, storage, model, dt, forces) where D
    state = storage.state
    conn_pos, conn_data = dispersion_cell_connections(model)
    update_total_darcy_velocity!(v, state.TotalVolumeFlux, state.ThermalDispersionGeometry,
        state.BulkVolume, conn_pos, conn_data, model.context, Val(D))
    return v
end

function update_total_darcy_velocity!(v, q, G, volumes, conn_pos, conn_data, context, dim::Val{D}) where D
    function F(cell)
        acc = zero(SVector{D, Float64})
        @inbounds for i in conn_pos[cell]:(conn_pos[cell + 1] - 1)
            conn = conn_data[i]
            face = conn.face
            C_l, C_r, = dispersion_face_geometry(G, face, dim)
            # Outward flux times vector from cell centroid to face centroid
            if conn.face_sign > 0
                acc += q[face]*C_l
            else
                acc -= q[face]*C_r
            end
        end
        V = value(volumes[cell])
        @inbounds for d in 1:D
            v[d, cell] = acc[d]/V
        end
    end
    Jutul.threaded_loop_minbatch(F, size(v, 2), context)
    return v
end

"""
    ThermalDispersionTransmissibility()

Half-face transmissibilities of the dispersion tensor (row 1: left cell, row 2:
right cell), evaluated at the start of each time step. The face transmissibility
is the harmonic average of the half-face values scaled by `ρ c_p` of each cell.
"""
struct ThermalDispersionTransmissibility <: VectorVariables end
Jutul.values_per_entity(model, ::ThermalDispersionTransmissibility) = 2
Jutul.associated_entity(::ThermalDispersionTransmissibility) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::ThermalDispersionTransmissibility, symb)
    return zeros(2, number_of_faces(data_domain))
end

function Jutul.update_parameter_before_step!(θ_d, ::ThermalDispersionTransmissibility, storage, model, dt, forces)
    state = storage.state
    if !haskey(state, :ThermalDispersivity)
        fill!(θ_d, 0.0)
        return θ_d
    end
    N = dispersion_neighborship(physical_representation(model.domain))
    v = state.TotalDarcyVelocity
    update_dispersion_transmissibility!(θ_d, N, state.ThermalDispersionGeometry, v,
        state.ThermalDispersivity, model.context, Val(size(v, 1)))
    return θ_d
end

function update_dispersion_transmissibility!(θ_d, N, G, v, α, context, dim::Val{D}) where D
    function F(face)
        @inbounds begin
            l, r = N[1, face], N[2, face]
            C_l, C_r, NA = dispersion_face_geometry(G, face, dim)
            θ_d[1, face] = dispersion_half_face_transmissibility(C_l, NA, v, α, l, dim)
            θ_d[2, face] = dispersion_half_face_transmissibility(C_r, -NA, v, α, r, dim)
        end
    end
    Jutul.threaded_loop_minbatch(F, size(θ_d, 2), context)
    return θ_d
end

"""
Half-face transmissibility `(D C)⋅N A / |C|²` for cell `cell` with
`D = α_T |v| I + (α_L - α_T) v v'/|v|`, `C` the vector from the cell centroid to
the face centroid and `NA` the area-weighted outward normal.
"""
@inline function dispersion_half_face_transmissibility(C, NA, v, α, cell, ::Val{D}) where D
    v_c = SVector{D, Float64}(ntuple(i -> @inbounds(v[i, cell]), Val(D)))
    vmag = norm(v_c)
    if vmag == 0.0
        return 0.0
    end
    @inbounds α_L = value(α[1, cell])
    @inbounds α_T = value(α[2, cell])
    DC_N = α_T*vmag*dot(C, NA) + (α_L - α_T)/vmag*dot(v_c, C)*dot(v_c, NA)
    return DC_N/dot(C, C)
end
