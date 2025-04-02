function get_multimodel_residuals(report)

    models = keys(report)

    residuals = Dict()

    for model in models
        rsd = get_model_residuals(report[model])
        residuals[model] = rsd
    end

    return residuals

end

function get_model_residuals(report)

    residuals = Dict()

    for eq_report in report

        equation = eq_report[:name]
        criterions = eq_report[:criterions]
        tolerances = eq_report[:tolerances]
        residual_norms = keys(criterions)

        equation_residuals = Dict()
    
        for res_norm in residual_norms

            rsd = criterions[res_norm].errors
            tol = tolerances[res_norm]
            nms = criterions[res_norm].names
            for (r, ϵ, α) in zip(rsd, tol, nms)
                α = process_name(α)
                if !haskey(equation_residuals, α)
                    equation_residuals[α] = Dict()
                end
                equation_residuals[α][res_norm] = r/ϵ
            end

        end

        residuals[equation] = equation_residuals
    
    end

    return residuals

end

function process_name(name)

    name = String(name)
    name = replace(name, " " => "_", "(" => "", ")" => "")
    name = Symbol(name)

end

function flatten_dict(input_dict::Dict, separator::String = ".", trail = [])
    values = []
    names = []

    for (key, value) in input_dict
        current_trail = vcat(trail, string(key))
        if value isa Dict
            sub_values, sub_names = flatten_dict(value, separator, current_trail)
            append!(values, sub_values)
            append!(names, sub_names)
        else
            push!(values, value)
            push!(names, join(current_trail, separator))
        end
    end

    return values, names
end