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
    cell,
    pressure = DEFAULT_MINIMUM_PRESSURE,
    temperature = 298.15;
    fractional_flow = nothing,
    density = nothing,
    trans_flow = 1e-12,
    trans_thermal = 1e-6
    )
    @assert isnothing(density) || density > 0.0 "Density, if provided, must be positive"
    if isnothing(fractional_flow)
        f = fractional_flow
    else
        f = Tuple(fractional_flow)
        @assert all(f .>= 0) "Fractional flow must be non-negative: $f"
        @assert sum(f) == 1.0 "Fractional flow for boundary condition in cell $cell must sum to 1."
    end
    @assert pressure >= DEFAULT_MINIMUM_PRESSURE
    @assert temperature >= 0.0
    return FlowBoundaryCondition(cell, pressure, temperature, trans_flow, trans_thermal, f, density)
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
    p = state.Pressure
    for bc in force
        c = bc.cell
        T_f = bc.trans_flow
        Δp = p[c] - bc.pressure
        q = T_f*Δp
        acc_i = view(acc, :, c)
        apply_flow_bc!(acc_i, q, bc, model, state, time)
    end
end

function apply_flow_bc!(acc, q, bc, model::SimulationModel{<:Any, T}, state, time) where T<:Union{ImmiscibleSystem, SinglePhaseSystem}
    mu = state.PhaseViscosities
    kr = state.RelativePermeabilities
    rho = state.PhaseMassDensities
    nph = length(acc)
    @assert size(kr, 1) == nph

    rho_inj = bc.density
    f_inj = bc.fractional_flow
    c = bc.cell
    if q > 0
        # Pressure inside is higher than outside, flow out from domain
        for ph in eachindex(acc)
            # Immiscible: Density * total flow rate * mobility for each phase
            q_i = q*rho[ph, c]*kr[ph, c]/mu[ph, c]
            acc[ph] += q_i
        end
    else
        # Injection of mass
        λ_t = 0.0
        for ph in eachindex(acc)
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
                acc[ph] += q*rho_inj*λ_t*F
            end
        else
            @assert length(f_inj) == nph
            for ph in 1:nph
                F = f_inj[ph]
                acc[ph] += q*rho_inj*λ_t*F
            end
        end
    end
end

function Jutul.vectorization_length(bc::FlowBoundaryCondition, variant)
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

function Jutul.vectorize_force!(v, bc::FlowBoundaryCondition, variant)
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

function Jutul.devectorize_force(bc::FlowBoundaryCondition, X, meta, variant)
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
