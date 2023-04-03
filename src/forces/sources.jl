"""
    SourceTerm(cell, value; fractional_flow = [1.0], type = MassSource)

Create source term in given `cell` with given total `value`. The optional
`fractional_flow` argument controls how this term is divided over components if
used for inflow.
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

