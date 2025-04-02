function compute_distance(report; 
    distance_function = r -> scaled_residual_norm(r),
    mapping = v -> maximum(v)
    )
    
    # Compute distance using distance function
    distance, names = distance_function(report)
    # Apply transformation
    distance = mapping(distance)

    return distance, names

end

function scaled_residual_norm(report)

    residuals = get_multimodel_residuals(report)
    values, names = flatten_dict(residuals)
    distance = max.(values .- 1.0, 0.0)

    return distance, names

end

function log_scaled_residual_norm(report)

    distance, names = scaled_residual_norm(report)
    distance = log10.(distance .+ 1.0)

    return distance, names

end

function nonconverged_equations(report)

    values, names = scaled_residual_norm(report)
    distance = (values .> 0.0) .+ 0.0

    return distance, names

end