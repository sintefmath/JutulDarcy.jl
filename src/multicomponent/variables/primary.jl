abstract type CompositionalFractions <: FractionVariables end

function values_per_entity(model, v::CompositionalFractions)
    sys = model.system
    nc = number_of_components(sys)
    if has_other_phase(sys)
        nval = nc - 1
    else
        nval = nc
    end
    return nval
end

function update_primary_variable!(state, p::CompositionalFractions, state_symbol, model, dx, w)
    s = state[state_symbol]
    Jutul.unit_sum_update!(s, p, model, dx, w)
end

struct OverallMoleFractions <: CompositionalFractions
    dz_max::Float64
end

"""
    OverallMoleFractions(;dz_max = 0.2)

Overall mole fractions definition for compositional. `dz_max` is the maximum
allowable change in any composition during a single Newton iteration.
"""
function OverallMoleFractions(;dz_max = 0.2)
    OverallMoleFractions(dz_max)
end
minimum_value(::OverallMoleFractions) = MultiComponentFlash.MINIMUM_COMPOSITION
absolute_increment_limit(z::OverallMoleFractions) = z.dz_max


function Jutul.increment_norm(dX, state, model, X, pvar::OverallMoleFractions)
    if haskey(state, :ImmiscibleSaturation)
        sw = state.ImmiscibleSaturation
    else
        sw = missing
    end
    T = eltype(dX)
    scale = @something Jutul.variable_scale(pvar) one(T)
    N = degrees_of_freedom_per_entity(model, pvar)
    sum_v = sum(abs, dX)
    max_v = maximum(abs, dX)
    sum_v_scaled, max_v_scaled = compositional_increment_scaled(
        dX, sw, model, Val(N))
    return (sum = scale*sum_v, sum_scaled = sum_v_scaled, max = scale*max_v, max_scaled = max_v_scaled)
end

function compositional_increment_scaled(dX, ::Missing, model, ::Val{N}) where N
    return sum(abs, dX), maximum(abs, dX)
end

function compositional_increment_scaled(dX, sw, model, ::Val{N}) where N
    M = global_map(model.domain)
    sw = Jutul.active_view(sw, M, for_variables = false)
    component_sums = ntuple(Val(N)) do component
        function scaled_increment(dx, saturation)
            return abs(dx)*(1.0 - value(saturation))
        end
        mapreduce(scaled_increment, +, view(dX, component, :), sw)
    end
    component_maxima = ntuple(Val(N)) do component
        function scaled_increment(dx, saturation)
            return abs(dx)*(1.0 - value(saturation))
        end
        mapreduce(scaled_increment, max, view(dX, component, :), sw)
    end
    return sum(component_sums), maximum(component_maxima)
end

"""
A single saturation variable that represents the "other" phase in a three phase
compositional system where two phases are predicted by an EoS
"""
Base.@kwdef struct ImmiscibleSaturation <: ScalarVariable
    ds_max::Float64 = 0.2
end

maximum_value(::ImmiscibleSaturation) = 1.0
minimum_value(::ImmiscibleSaturation) = 0.0
absolute_increment_limit(s::ImmiscibleSaturation) = s.ds_max
