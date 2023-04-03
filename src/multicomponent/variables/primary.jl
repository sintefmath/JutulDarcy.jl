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

export OverallMoleFractions
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


"""
A single saturation that represents the "other" phase in a
three phase compositional system where two phases are predicted by an EoS
"""
Base.@kwdef struct ImmiscibleSaturation <: ScalarVariable
    ds_max::Float64 = 0.2
end

maximum_value(::ImmiscibleSaturation) = 1.0 - MINIMUM_COMPOSITIONAL_SATURATION
minimum_value(::ImmiscibleSaturation) = 0.0
absolute_increment_limit(s::ImmiscibleSaturation) = s.ds_max
