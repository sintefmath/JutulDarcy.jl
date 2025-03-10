function Jutul.vectorization_length(v, model::FacilityModel, controls_or_limits::Dict, name, variant)
    error("! $name")
    error("Not implemented for $(typeof(forces))")
end

function Jutul.vectorize_force!(v, model::FacilityModel, controls_or_limits::Dict, name, variant)
    error("!")
end
