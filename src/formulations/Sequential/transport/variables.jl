struct TotalSaturation <: Jutul.ScalarVariable
    max_abs::Union{Float64, Nothing}
    max_rel::Union{Float64, Nothing}
    minimum_value::Float64
    maximum_value::Float64
end

function TotalSaturation(;
        max_abs = nothing,
        max_rel = 0.2,
        minimum_value = 1e-3,
        maximum_value = Inf
    )
    if !isnothing(max_abs)
        max_abs > 0.0 || throw(ArgumentError("Maximum absolute total saturation change was $max_abs, must be positive"))
    end
    if !isnothing(max_rel)
        max_rel > 0.0 || throw(ArgumentError("Maximum relative total saturation change was $max_rel, must be positive"))
    end
    minimum_value < maximum_value || throw(ArgumentError("Maximum pressure $maximum_value must be larger than minimum pressure $minimum_pressure"))
    TotalSaturation(max_abs, max_rel, minimum_value, maximum_value)
end

Jutul.minimum_value(st::TotalSaturation) = st.minimum_value
Jutul.maximum_value(st::TotalSaturation) = st.maximum_value
Jutul.absolute_increment_limit(st::TotalSaturation) = st.max_abs
Jutul.relative_increment_limit(st::TotalSaturation) = st.max_rel

struct TotalVolumetricFlux <: Jutul.ScalarVariable end
Jutul.associated_entity(::TotalVolumetricFlux) = Faces()

struct PerforationTotalVolumetricFlux <: Jutul.ScalarVariable end
Jutul.associated_entity(::PerforationTotalVolumetricFlux) = JutulDarcy.Perforations()


struct TotalSaturationCorrectedScalarVariable <: Jutul.ScalarVariable
    uncorrected_label::Symbol
end

struct TotalSaturationCorrectedVariable <: Jutul.VectorVariables
    uncorrected_label::Symbol
    values_per_entity::Int
end

Jutul.values_per_entity(model, ts::TotalSaturationCorrectedVariable) = ts.values_per_entity

function Jutul.get_dependencies(ts::Union{TotalSaturationCorrectedVariable, TotalSaturationCorrectedScalarVariable}, model::TransportModel)
    return (:TotalSaturation, ts.uncorrected_label)
end

function Jutul.update_secondary_variable!(v, ts::Union{TotalSaturationCorrectedVariable, TotalSaturationCorrectedScalarVariable}, model, state, ix = entity_eachindex(v))
    # Get the uncorrected variable
    UncorrectedVariable = state[ts.uncorrected_label]
    # Get the total saturation
    TotalSaturation = state.TotalSaturation
    # Update the total saturation corrected variable
    total_sat_correction(v, ts, model, UncorrectedVariable, TotalSaturation, ix)
end

function total_sat_correction(v, ts::TotalSaturationCorrectedVariable, model::TransportModel, UncorrectedVariable, TotalSaturation, ix)
    for i in ix
        sT = TotalSaturation[i]
        for j in axes(v, 1)
            v[j, i] = sT*UncorrectedVariable[j, i]
        end
    end
end

function total_sat_correction(v, ts::TotalSaturationCorrectedScalarVariable, model::TransportModel, UncorrectedVariable, TotalSaturation, ix)
    for i in ix
        sT = TotalSaturation[i]
        v[i] = sT*UncorrectedVariable[i]
    end
end
