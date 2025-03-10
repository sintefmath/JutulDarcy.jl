"""
    SourceTerm(cell, value; fractional_flow = [1.0], type = MassSource)

Create source term in given `cell` with given total `value`.

The optional `fractional_flow` argument controls how this term is divided over
components if used for inflow and should contain one entry per component in the
system: (`number_of_components(system)`). `fractional_flow` should sum up to
1.0. The `type` argument should be an instance of the `FlowSourceType` enum,
with interpretations as follows:

- `MassSource`: Source is directly interpreted as component masses.
- `StandardVolumeSource`: Source is volume at standard/surface conditions.
   References densities are used to convert into mass sources.
- `VolumeSource`: Source is volume at in-situ / reservoir conditions.

"""
function SourceTerm(cell, value; fractional_flow = [1.0], type = MassSource)
    @assert sum(fractional_flow) == 1.0 "Fractional flow for source term in cell $cell must sum to 1."
    f = Tuple(fractional_flow)
    return SourceTerm(cell, value, f, type)
end

function cell(s::SourceTerm{I, T}) where {I, T} 
    return s.cell::I
end

function Jutul.subforce(s::AbstractVector{S}, model) where S<:SourceTerm
    # Just to be safe
    s = deepcopy(s)
    m = global_map(model.domain)

    n = length(s)
    keep = repeat([false], n)
    for (i, src) in enumerate(s)
        # Cell must be in local domain, and not on boundary
        if !Jutul.global_cell_inside_domain(src.cell, m)
            continue
        end
        c_l = Jutul.local_cell(src.cell, m)
        c_i = Jutul.interior_cell(c_l, m)
        inner = !isnothing(c_i)
        if !inner
            continue
        end
        keep[i] = true
        s[i] = SourceTerm(c_i, src.value, fractional_flow = src.fractional_flow)
    end
    return s[keep]
end

function Jutul.vectorization_length(src::SourceTerm, model, name, variant)
    n = 1
    if variant == :all
        f = src.fractional_flow
        if !isnothing(f)
            n += length(f)
        end
    elseif variant == :control
        # Do nothing
    else
        error("Variant $variant not supported")
    end
    return n
end

function Jutul.vectorize_force!(v, model::SimulationModel, src::SourceTerm, name, variant)
    v[1] = src.value
    names = [:value]
    if variant == :all
        offset = 1
        f = src.fractional_flow
        if !isnothing(f)
            for (i, f_i) in enumerate(f)
                offset += 1
                v[offset] = f_i
                push!(names, Symbol("fractional_flow$i"))
            end
        end
    elseif variant == :control
        # Do nothing
    else
        error("Variant $variant not supported")
    end
    return (names = names, )
end

function Jutul.devectorize_force(src::SourceTerm, model::SimulationModel, X, meta, name, variant)
    val = X[1]
    f = src.fractional_flow
    if variant == :all
        if !isnothing(f)
            f = tuple(X[2:end]...)
        end
    elseif variant == :control
        # Do nothing
    else
        error("Variant $variant not supported")
    end
    return SourceTerm(src.cell, val, f, src.type)
end
