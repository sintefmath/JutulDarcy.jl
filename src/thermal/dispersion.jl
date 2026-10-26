"""
    add_thermal_dispersion!(domain::DataDomain, α)

Add thermal dispersivity to a reservoir domain. The dispersivity `α` (in meters)
can be given as

- a number, used as both longitudinal and transverse dispersivity in all cells,
- a tuple or vector `(α_L, α_T)` with longitudinal and transverse dispersivity
  for all cells,
- a vector with one value per cell, used for both components,
- a matrix of size `2 × num_cells` with `α_L` in the first row and `α_T` in the
  second row.

The dispersive heat flux is `-Σ_α ρ_α c_α D_α ∇T` with
`D_α = α_T |v_α| I + (α_L - α_T) v_α v_α'/|v_α|`, where `v_α` is the Darcy
velocity of phase `α` and `c_α` the specific heat capacity of the phase (see
`FluidHeatCapacity`). The dispersion enters the energy equation as an
additional face transmissibility. Note that:

- `D_α` is evaluated explicitly from the fluxes at the start of each time step.
  The dependence of the dispersion on the flow is therefore not part of the
  Jacobian, and sensitivities with respect to the dispersivity (or other
  parameters through the velocity) are not computed by the adjoint method.
- Thermal dispersion requires a two-point flux discretization, where only the
  projection of `D_α` onto the face normal is used. The cross terms of the
  tensor are therefore not represented, and the discretization is not
  consistent on grids that are not K-orthogonal for the dispersion tensor.
"""
function add_thermal_dispersion!(domain::DataDomain, α)
    domain[:thermal_dispersivity, Cells()] = thermal_dispersivity_values(α, number_of_cells(domain))
    return domain
end

"""
    thermal_dispersivity_values(α, nc)

Convert the input to [`add_thermal_dispersion!`](@ref) to a `2 × nc` matrix with
longitudinal (first row) and transverse (second row) dispersivity.
"""
function thermal_dispersivity_values(α, nc)
    if α isa Number
        D = fill(Float64(α), 2, nc)
    elseif α isa Tuple
        length(α) == 2 || throw(ArgumentError("Expected tuple (α_L, α_T) for thermal dispersivity, got $α"))
        D = repeat([Float64(α[1]), Float64(α[2])], 1, nc)
    elseif α isa AbstractVector
        n = length(α)
        if n == 2
            D = repeat(Float64.(α), 1, nc)
        elseif n == nc
            D = repeat(Float64.(α)', 2, 1)
        else
            throw(ArgumentError("Expected thermal dispersivity vector of length 2 or num_cells ($nc), got length $n"))
        end
    elseif α isa AbstractMatrix
        size(α) == (2, nc) || throw(ArgumentError("Expected thermal dispersivity matrix of size (2, $nc), got size $(size(α))"))
        D = Float64.(α)
    else
        throw(ArgumentError("Unsupported type $(typeof(α)) for thermal dispersivity"))
    end
    all(x -> isfinite(x) && x >= 0, D) || throw(ArgumentError("Thermal dispersivity must be finite and non-negative"))
    return D
end

"""
    ThermalDispersivity()

Cell-wise thermal dispersivity (m) with two components per cell: row 1 is
longitudinal and row 2 is transverse. Initialized from `:thermal_dispersivity`
in the domain, see [`add_thermal_dispersion!`](@ref).
"""
struct ThermalDispersivity <: VectorVariables end
Jutul.values_per_entity(model, ::ThermalDispersivity) = 2
Jutul.minimum_value(::ThermalDispersivity) = 0.0
Jutul.associated_entity(::ThermalDispersivity) = Cells()

function Jutul.default_parameter_values(data_domain, model, param::ThermalDispersivity, symb)
    haskey(data_domain, :thermal_dispersivity) || error(":thermal_dispersivity must be present in DataDomain to initialize parameter $symb")
    return thermal_dispersivity_values(data_domain[:thermal_dispersivity], number_of_cells(data_domain))
end

# Dispersion is computed in three stages before each time step, all in the
# update of ThermalDispersionTransmissibility:
#   1. PhaseVolumeFluxes: Darcy volume flux of each phase over each face
#   2. PhaseDarcyVelocities: cell velocities reconstructed from the face fluxes
#   3. ThermalDispersionTransmissibility: half-face transmissibilities of D
# All topology is taken from the model's own discretization and all geometry
# from the face parameter ThermalDispersionGeometry, so that the computations
# only use data that is local to the model.

function add_thermal_dispersion_to_model!(model)
    model_uses_two_point_flow(model) || throw(ArgumentError("Thermal dispersion requires a two-point flux discretization."))
    dim = size(model.data_domain[:cell_centroids], 1)
    set_parameters!(model,
        ThermalDispersivity = ThermalDispersivity(),
        ThermalDispersionGeometry = ThermalDispersionGeometry(dim),
        PhaseVolumeFluxes = PhaseVolumeFluxes(),
        PhaseDarcyVelocities = PhaseDarcyVelocities(dim),
        ThermalDispersionTransmissibility = ThermalDispersionTransmissibility(),
    )
    set_secondary_variables!(model,
        FluidHeatCapacity = FluidHeatCapacity()
    )
    return model
end

"""
    FluidHeatCapacity()

Specific heat capacity at constant pressure of each fluid phase (J/kg/K), used
for thermal dispersion. The default takes the value from
`ComponentHeatCapacity`, which requires one component per phase or a
compositional model. Models where the enthalpy is not given by
`ComponentHeatCapacity`, such as tabulated enthalpies, should set this
variable explicitly, e.g. as a
`PressureTemperatureDependentVariable` from a heat capacity table.
"""
struct FluidHeatCapacity <: PhaseVariables end
Jutul.minimum_value(::FluidHeatCapacity) = 0.0

@jutul_secondary function update_fluid_heat_capacity!(c_p, var::FluidHeatCapacity, model, ComponentHeatCapacity, ix)
    C = ComponentHeatCapacity
    @assert size(c_p) == size(C) "FluidHeatCapacity from ComponentHeatCapacity assumes one component per phase. Set FluidHeatCapacity explicitly for this model."
    for i in ix
        for ph in axes(c_p, 1)
            c_p[ph, i] = C[ph, i]
        end
    end
    return c_p
end

@jutul_secondary function update_fluid_heat_capacity!(c_p, var::FluidHeatCapacity, model::CompositionalModel, ComponentHeatCapacity, LiquidMassFractions, VaporMassFractions, ix)
    fsys = model.system
    C = ComponentHeatCapacity
    X = LiquidMassFractions
    Y = VaporMassFractions
    if has_other_phase(fsys)
        a, l, v = phase_indices(fsys)
        offset = 1
    else
        l, v = phase_indices(fsys)
        offset = 0
    end
    for i in ix
        if has_other_phase(fsys)
            c_p[a, i] = C[1, i]
        end
        C_l = 0.0
        C_v = 0.0
        for c in axes(X, 1)
            C_i = C[c+offset, i]
            C_l += C_i*X[c, i]
            C_v += C_i*Y[c, i]
        end
        c_p[l, i] = C_l
        c_p[v, i] = C_v
    end
    return c_p
end

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
    PhaseVolumeFluxes()

Darcy volume flux of each phase over each face (positive from the left to the
right cell of the face), evaluated at the start of each time step when
`ThermalDispersionTransmissibility` is updated.
"""
struct PhaseVolumeFluxes <: VectorVariables end
Jutul.values_per_entity(model, ::PhaseVolumeFluxes) = number_of_phases(model.system)
Jutul.associated_entity(::PhaseVolumeFluxes) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::PhaseVolumeFluxes, symb)
    return zeros(number_of_phases(model.system), number_of_faces(data_domain))
end

function update_phase_volume_fluxes!(q, N, state, model)
    function F(face)
        @inbounds l, r = N[1, face], N[2, face]
        q_face = darcy_phase_volume_fluxes(face, state, model, Jutul.DefaultFlux(), TPFA(l, r, 1), SPU(l, r))
        for ph in axes(q, 1)
            @inbounds q[ph, face] = value(q_face[ph])
        end
    end
    Jutul.threaded_loop_minbatch(F, size(q, 2), model.context)
    return q
end

"""
    PhaseDarcyVelocities(dim)

Cell-wise Darcy velocity of each phase reconstructed from the face fluxes as
`v_c = Σ_f q_f (x_f - x_c) / V_c`, evaluated at the start of each time step when
`ThermalDispersionTransmissibility` is updated.
Rows `(α-1)*dim+1:α*dim` hold the velocity of phase `α`. Sources (wells and
boundary conditions) are located at the cell centroid and do not contribute.
The normal component is replaced by the face flux for each half-face when the
dispersion tensor is evaluated (see `ThermalDispersionTransmissibility`).
"""
struct PhaseDarcyVelocities{D} <: VectorVariables end
PhaseDarcyVelocities(dim::Int) = PhaseDarcyVelocities{dim}()
Jutul.values_per_entity(model, ::PhaseDarcyVelocities{D}) where D = D*number_of_phases(model.system)
Jutul.associated_entity(::PhaseDarcyVelocities) = Cells()

function Jutul.default_parameter_values(data_domain, model, param::PhaseDarcyVelocities{D}, symb) where D
    return zeros(D*number_of_phases(model.system), number_of_cells(data_domain))
end

function update_phase_darcy_velocities!(v, q, G, volumes, conn_pos, conn_data, context, dim::Val{D}) where D
    nph = size(q, 1)
    function F(cell)
        V = value(volumes[cell])
        for ph in 1:nph
            acc = zero(SVector{D, Float64})
            @inbounds for i in conn_pos[cell]:(conn_pos[cell + 1] - 1)
                conn = conn_data[i]
                face = conn.face
                C_l, C_r, = dispersion_face_geometry(G, face, dim)
                # Outward flux times vector from cell centroid to face centroid
                if conn.face_sign > 0
                    acc += q[ph, face]*C_l
                else
                    acc -= q[ph, face]*C_r
                end
            end
            offset = (ph - 1)*D
            @inbounds for d in 1:D
                v[offset + d, cell] = acc[d]/V
            end
        end
    end
    Jutul.threaded_loop_minbatch(F, size(v, 2), context)
    return v
end

"""
    ThermalDispersionTransmissibility()

Half-face transmissibilities of the dispersion tensor of each phase (row
`2α-1`: left cell, row `2α`: right cell), evaluated at the start of each time
step. The face transmissibility is the harmonic average of the half-face values
scaled by `ρ_α c_α` of each cell. The transmissibility is zero for faces where
the static fluid thermal transmissibility of the phase is zero (e.g. faces
blocked by fractures).
"""
struct ThermalDispersionTransmissibility <: VectorVariables end
Jutul.values_per_entity(model, ::ThermalDispersionTransmissibility) = 2*number_of_phases(model.system)
Jutul.associated_entity(::ThermalDispersionTransmissibility) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::ThermalDispersionTransmissibility, symb)
    return zeros(2*number_of_phases(model.system), number_of_faces(data_domain))
end

function Jutul.update_parameter_before_step!(θ_d, ::ThermalDispersionTransmissibility, storage, model, dt, forces)
    state = storage.state
    N = dispersion_neighborship(physical_representation(model.domain))
    G = state.ThermalDispersionGeometry
    dim = Val(size(G, 1) ÷ 3)
    # Fluxes and velocities are updated here, so that the result does not
    # depend on the order of the parameters.
    q = state.PhaseVolumeFluxes
    v = state.PhaseDarcyVelocities
    update_phase_volume_fluxes!(q, N, state, model)
    conn_pos, conn_data = dispersion_cell_connections(model)
    update_phase_darcy_velocities!(v, q, G, state.BulkVolume, conn_pos, conn_data, model.context, dim)
    update_dispersion_transmissibility!(θ_d, N, G, q, v, state.ThermalDispersivity,
        state.FluidThermalConductivities, model.context, dim)
    return θ_d
end

function update_dispersion_transmissibility!(θ_d, N, G, q, v, α, θ_f, context, dim::Val{D}) where D
    nph = size(q, 1)
    function F(face)
        @inbounds begin
            l, r = N[1, face], N[2, face]
            C_l, C_r, NA = dispersion_face_geometry(G, face, dim)
            for ph in 1:nph
                if θ_f[ph, face] > 0
                    q_ph = q[ph, face]
                    v_l = dispersion_cell_velocity(v, ph, l, dim)
                    v_r = dispersion_cell_velocity(v, ph, r, dim)
                    θ_l = dispersion_half_face_transmissibility(C_l, NA, q_ph, v_l, α, l)
                    θ_r = dispersion_half_face_transmissibility(C_r, -NA, -q_ph, v_r, α, r)
                else
                    θ_l = θ_r = 0.0
                end
                θ_d[2ph - 1, face] = θ_l
                θ_d[2ph, face] = θ_r
            end
        end
    end
    Jutul.threaded_loop_minbatch(F, size(θ_d, 2), context)
    return θ_d
end

@inline function dispersion_cell_velocity(v, ph, cell, ::Val{D}) where D
    offset = (ph - 1)*D
    return SVector{D, Float64}(ntuple(i -> @inbounds(v[offset + i, cell]), Val(D)))
end

"""
Half-face transmissibility `max((D C)⋅N A / |C|², 0)` for cell `cell` with
`D = α_T |v| I + (α_L - α_T) v v'/|v|`, `C` the vector from the cell centroid to
the face centroid and `NA` the area-weighted outward normal. The velocity `v` is
the cell velocity `v_c` with its normal component replaced by the outward face
flux `q_out` divided by the face area, so that the flow over each face (e.g.
radial flow around a well) is represented even if the cell average vanishes.
Negative values, which can occur for anisotropic dispersion on grids that are
not K-orthogonal, are set to zero.
"""
@inline function dispersion_half_face_transmissibility(C, NA, q_out, v_c, α, cell)
    v = v_c + (q_out - dot(NA, v_c))/dot(NA, NA)*NA
    vmag = norm(v)
    if vmag == 0.0
        return 0.0
    end
    @inbounds α_L = value(α[1, cell])
    @inbounds α_T = value(α[2, cell])
    DC_N = α_T*vmag*dot(C, NA) + (α_L - α_T)/vmag*dot(v, C)*dot(v, NA)
    return max(DC_N/dot(C, C), 0.0)
end

"""
    thermal_dispersion_face_transmissibility(state, face, tpfa, nph)

Thermal dispersion transmissibility of a face, summed over the phases. For each
phase, the half-face transmissibilities are scaled by `ρ_α c_α` of the cell and
combined with a harmonic average.
"""
@inline function thermal_dispersion_face_transmissibility(state, face, tpfa::TPFA, nph)
    θ_d = state.ThermalDispersionTransmissibility
    ρ = state.PhaseMassDensities
    c_p = state.FluidHeatCapacity
    l, r = tpfa_face_cells(tpfa)
    θ = zero(eltype(ρ))*zero(eltype(c_p))
    for ph in 1:nph
        @inbounds θ_l = θ_d[2ph - 1, face]
        @inbounds θ_r = θ_d[2ph, face]
        if θ_l > 0 && θ_r > 0
            a_l = @inbounds θ_l*ρ[ph, l]*c_p[ph, l]
            a_r = @inbounds θ_r*ρ[ph, r]*c_p[ph, r]
            θ += a_l*a_r/(a_l + a_r)
        end
    end
    return θ
end
