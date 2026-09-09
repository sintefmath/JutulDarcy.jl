# Structural adaptation for the compact topology used by the single-model
# reservoir discretization. More complex property types can add equivalent
# Adapt rules without changing the execution context.
function Adapt.adapt_structure(to, g::MinimalTPFATopology)
    return MinimalTPFATopology(g.nc, Adapt.adapt(to, g.neighborship))
end
