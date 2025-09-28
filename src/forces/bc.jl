"""
    FlowBoundaryCondition(
    cell,
    pressure = DEFAULT_MINIMUM_PRESSURE,
    temperature = 298.15;
    fractional_flow = nothing,
    density = nothing,
    trans_flow = 1e-12,
    trans_thermal = 1e-6
    )

Dirchlet boundary condition for constant values (pressure/temperature) at some inflow boundary
"""
function FlowBoundaryCondition(
        cell::Int,
        pressure = DEFAULT_MINIMUM_PRESSURE,
        temperature = 298.15;
        fractional_flow = nothing,
        density = nothing,
        trans_flow = 1e-12,
        trans_thermal = 1e-6
    )
    pressure >= DEFAULT_MINIMUM_PRESSURE || throw(ArgumentError("Pressure must be at least $DEFAULT_MINIMUM_PRESSURE"))
    temperature >= 0.0 || throw(ArgumentError("Temperature must be at least 0.0 K"))

    isnothing(density) || density > 0.0 || error("Density, if provided, must be positive")
    if isnothing(fractional_flow)
        f = fractional_flow
    else
        T = promote_type(eltype(fractional_flow), typeof(pressure))
        fractional_flow = convert.(T, fractional_flow)
        pressure = convert(T, pressure)
        f = Tuple(fractional_flow)
        all(f .>= 0) || error("Fractional flow must be non-negative: $f")
        sum(f) == 1.0 || error("Fractional flow for boundary condition in cell $cell must sum to 1.")
    end
    pressure, temperature, trans_flow, trans_thermal = promote(pressure, temperature, trans_flow, trans_thermal)
    return FlowBoundaryCondition(cell, pressure, temperature, trans_flow, trans_thermal, f, density)
end

function FlowBoundaryCondition(
        domain::DataDomain,
        cell::Int,
        pressure = DEFAULT_MINIMUM_PRESSURE,
        temperature = 298.15;
        dir = :x,
        kwarg...
    )
    G = physical_representation(domain)
    D = Jutul.dim(G)
    dir in (:x, :y, :z) || dir in 1:3 || throw(ArgumentError("Direction argument `dir` must be either :x, :y, :z or 1, 2, 3"))
    G = physical_representation(domain)
    if dir isa Symbol
        dir = findfirst(isequal(dir), (:x, :y, :z))
    end
    dist = cell_dims(G, cell)[dir]
    # Approximate area since we don't know the face.
    A = 1.0
    for i in 1:D
        if i != dir
            A *= cell_dims(G, cell)[i]
        end
    end
    K = domain[:permeability]
    cond = domain[:rock_thermal_conductivity]

    function local_trans(perm_or_c)
        if perm_or_c isa Vector
            ki = perm_or_c[cell]
        else
            ki = perm_or_c[:, cell]
        end
        # Take the diagonal
        ki = Jutul.expand_perm(ki, D)[dir, dir]
        # Distance to boundary is half the cell width
        return A*ki/(dist/2.0)
    end

    T_flow = local_trans(K)
    T_heat = local_trans(cond)
    return FlowBoundaryCondition(cell, pressure, temperature;
        trans_flow = T_flow,
        trans_thermal = T_heat,
        kwarg...
    )
end

export flow_boundary_condition
"""
    flow_boundary_condition(cells, domain, pressures, temperatures = 298.15; kwarg...)

Add flow boundary conditions to a vector of `cells` for a given `domain` coming
from `reservoir_domain`. The input arguments `pressures` and `temperatures` can
either be scalars or one value per cell. Other keyword arguments are passed onto
the `FlowBoundaryCondition` constructor.

The output of this function is a `Vector` of boundary conditions that can be
passed on the form `forces = setup_reservoir_forces(model, bc = bc)`.
"""
function flow_boundary_condition(cells, domain, pressures, temperatures = 298.15; fractional_flow = nothing, kwarg...)
    if fractional_flow isa Vector
        fractional_flow = tuple(fractional_flow...)
    end
    bc = []
    flow_boundary_condition!(bc, domain, cells, pressures, temperatures; fractional_flow = fractional_flow, kwarg...)
    return [i for i in bc]
end

function flow_boundary_condition!(bc, domain, cells, pressures, temperatures = 298.15; kwarg...)
    n = length(cells)
    if temperatures isa Real
        temperatures = fill(temperatures, n)
    end
    if pressures isa Real
        pressures = fill(pressures, n)
    end
    length(temperatures) == n || throw(ArgumentError("Mismatch in length of cells and temperatures arrays"))
    length(pressures) == n || throw(ArgumentError("Mismatch in length of cells and pressures arrays"))

    for (cell, pressure, temperature) in zip(cells, pressures, temperatures)
        bc_c = FlowBoundaryCondition(domain, cell, pressure, temperature; kwarg...)
        push!(bc, bc_c)
    end

    return bc
end

function Jutul.subforce(s::AbstractVector{S}, model) where S<:FlowBoundaryCondition
    s = deepcopy(s)
    m = global_map(model.domain)
    n = length(s)
    keep = repeat([false], n)
    for (i, bc) in enumerate(s)
        # Cell must be in local domain, and not on boundary
        if !Jutul.global_cell_inside_domain(bc.cell, m)
            continue
        end
        c_l = Jutul.local_cell(bc.cell, m)
        c_i = Jutul.interior_cell(c_l, m)
        inner = !isnothing(c_i)
        if !inner
            continue
        end
        keep[i] = true
        s[i] = FlowBoundaryCondition(
            c_i,
            bc.pressure,
            bc.temperature,
            bc.trans_flow,
            bc.trans_thermal,
            bc.fractional_flow,
            bc.density
        )
    end
    return s[keep]
end

function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel{D, S}, eq::ConservationLaw{:TotalMasses}, eq_s, force::V, time) where {V <: AbstractVector{<:FlowBoundaryCondition}, D, S<:MultiPhaseSystem}
    state = storage.state
    nph = number_of_phases(reservoir_model(model).system)
    for bc in force
        c = bc.cell
        acc_i = view(acc, :, c)
        q = compute_bc_mass_fluxes(bc, global_map(model), state, nph)
        apply_flow_bc!(acc_i, q, bc, model, state, time)
    end
end

function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel{D, S}, eq::ConservationLaw{:TotalThermalEnergy}, eq_s, force::V, time) where {V <: AbstractVector{<:FlowBoundaryCondition}, D, S<:MultiPhaseSystem}
    state = storage.state
    nph = number_of_phases(reservoir_model(model).system)
    for bc in force
        c = bc.cell
        acc_i = view(acc, :, c)
        qh_adv, qh_cond = compute_bc_heat_fluxes(bc, global_map(model), state, nph)
        apply_flow_bc!(acc_i, qh_adv + qh_cond, bc, model, state, time)
    end
end

function compute_bc_mass_fluxes(bc, gmap, state, nph)
    # Get reservoir properties
    p   = state.Pressure
    mu  = state.PhaseViscosities
    kr  = state.RelativePermeabilities
    rho = state.PhaseMassDensities
    @assert size(kr, 1) == nph
    # Get boundary properties
    c       = Jutul.full_cell(bc.cell, gmap)
    T_f     = bc.trans_flow
    rho_inj = bc.density
    f_inj   = bc.fractional_flow
    # Compute total mass flux
    Δp = p[c] - bc.pressure
    q_tot = T_f*Δp

    num_t = Base.promote_type(typeof(q_tot), eltype(kr), eltype(mu), eltype(rho), typeof(q_tot))
    isbits_out = isbitstype(num_t)
    if isbits_out
        V_t = MVector{nph, num_t}
    else
        V_t = SizedVector{nph, num_t}
    end
    q = zeros(V_t)
    if q_tot > 0
        # Pressure inside is higher than outside, flow out from domain
        for ph in 1:nph
            # Immiscible: Density * total flow rate * mobility for each phase
            q[ph] = q_tot*rho[ph, c]*kr[ph, c]/mu[ph, c]
        end
    else
        # Injection of mass
        λ_t = 0.0
        for ph in 1:nph
            λ_t += kr[ph, c]/mu[ph, c]
        end
        if isnothing(rho_inj)
            # Density not provided, take saturation average from what we have in
            # the inside of the domain
            rho_inj = 0.0
            for ph in 1:nph
                rho_inj += state.Saturations[ph, c]*rho[ph, c]
            end
        end
        if isnothing(f_inj)
            # Fractional flow not provided. We match the mass fraction we
            # observe on the inside.
            total = 0.0
            for ph in 1:nph
                total += state.TotalMasses[ph, c]
            end
            for ph in 1:nph
                F = state.TotalMasses[ph, c]/total
                q[ph] = q_tot*rho_inj*λ_t*F
            end
        else
            @assert length(f_inj) == nph
            for ph in 1:nph
                F = f_inj[ph]
                q[ph] = q_tot*rho_inj*λ_t*F
            end
        end
    end
    if isbits_out
        out = SVector{nph, num_t}(q)
    else
        out = q
    end
    return q
end

function compute_bc_heat_fluxes(bc, gmap, state, nph)
    q = compute_bc_mass_fluxes(bc, gmap, state, nph)
    c = Jutul.full_cell(bc.cell, gmap)

    # Get reservoir properties
    p = state.Pressure
    T = state.Temperature
    s = state.Saturations
    h = state.FluidEnthalpy
    u = state.FluidInternalEnergy
    rho = state.PhaseMassDensities

    # Get boundary properties
    T_h    = bc.trans_thermal
    p_bc   = bc.pressure
    T_bc   = bc.temperature
    rho_bc = bc.density

    qh_advective = 0
    for ph in 1:nph
        if q[ph] > 0
            # Flow out from domain
            h_ph = h[ph,c]
        else
            Cp = u[ph,c]/T[c]
            if isnothing(rho_bc)
                # Density not provided, take saturation average from what we
                # have in the inside of the domain
                rho_bc = 0.0
                for ph in 1:nph
                    rho_bc += state.Saturations[ph,c]*rho[ph,c]
                end
            end
            h_ph = Cp*T_bc + p_bc/rho_bc
        end
        qh_advective += h_ph*q[ph]
    end

    ΔT = T[c] - T_bc
    qh_conductive = T_h*ΔT

    return qh_advective, qh_conductive
end

function apply_flow_bc!(acc, q, bc, model::SimulationModel{<:Any, T}, state, time) where T<:Union{ImmiscibleSystem, SinglePhaseSystem}

    for ph in eachindex(acc)
        acc[ph] += q[ph]
    end

end

function Jutul.vectorization_length(bc::FlowBoundaryCondition, model, name, variant)
    if variant == :all
        n = 4 # pressure, temperature, T_flow, T_thermal
        f = bc.fractional_flow
        if !isnothing(f)
            n += length(f)
        end
        if !isnothing(bc.density)
            n += 1
        end
        return n
    elseif variant == :control
        return 1
    else
        error("Variant $variant not supported")
    end
end

function Jutul.vectorize_force!(v, model::SimulationModel, bc::FlowBoundaryCondition, name, variant)
    names = []
    if variant == :all
        v[1] = bc.pressure
        push!(names, :pressure)
        v[2] = bc.temperature
        push!(names, :temperature)
        v[3] = bc.trans_flow
        push!(names, :trans_flow)
        v[4] = bc.trans_thermal
        push!(names, :trans_thermal)
        offset = 4
        f = bc.fractional_flow
        if !isnothing(f)
            for (i, f_i) in enumerate(f)
                offset += 1
                v[offset] = f_i
                push!(names, Symbol("fractional_flow$i"))
            end
        end
        if !isnothing(bc.density)
            offset += 1
            v[offset] = bc.density
            push!(names, :density)
        end
    elseif variant == :control
        v[1] = bc.pressure
        push!(names, :pressure)
    else
        error("Variant $variant not supported")
    end
    return (names = names, )
end

function Jutul.devectorize_force(bc::FlowBoundaryCondition, model::SimulationModel, X, meta, name, variant)
    p = X[1]
    T = bc.temperature
    trans_flow = bc.trans_flow
    trans_thermal = bc.trans_thermal
    f = bc.fractional_flow
    ρ = bc.density
    if variant == :all
        T = X[2]
        trans_flow = X[3]
        trans_thermal = X[4]

        offset = 4
        if !isnothing(f)
            f = X[(offset+1):(offset+length(f))]
        end
        if !isnothing(ρ)
            ρ = X[end]
        end
    elseif variant == :control
        # DO nothing
    else
        error("Variant $variant not supported")
    end
    return FlowBoundaryCondition(
        bc.cell,
        p,
        T,
        trans_flow,
        trans_thermal,
        f,
        ρ
    )
end
