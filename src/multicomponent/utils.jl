
get_components(sys::MultiComponentSystem) = sys.components
number_of_components(sys::MultiComponentSystem) = length(get_components(sys))

liquid_phase_index(sys::MultiPhaseCompositionalSystemLV) = phase_index(sys, LiquidPhase())
vapor_phase_index(sys::MultiPhaseCompositionalSystemLV) = phase_index(sys, VaporPhase())
other_phase_index(sys::MultiPhaseCompositionalSystemLV{E, T, O}) where {E, T, O} = phase_index(sys, O())

function number_of_components(sys::MultiPhaseCompositionalSystemLV{E, T, O, G, N}) where {E, T, O, G, N}
    return N + has_other_phase(sys)
end

phase_index(sys, phase) = only(findfirst(isequal(phase), sys.phases))
has_other_phase(sys) = number_of_phases(sys) > 2
has_other_phase(sys::CompositeSystem) = has_other_phase(flow_system(sys))
has_other_phase(sys::MultiPhaseCompositionalSystemLV{E, T, O}) where {E, T, O<:Nothing} = false

phase_indices(sys::MultiPhaseCompositionalSystemLV{E, T, O}) where {E, T, O<:Nothing} = (liquid_phase_index(sys), vapor_phase_index(sys))
phase_indices(sys::MultiComponentSystem) = (other_phase_index(sys), liquid_phase_index(sys), vapor_phase_index(sys))

export KValueWrapper
struct KValueWrapper{T, D}
    K::T
end

"""
    KValueWrapper(K; dependence::Symbol = :pT)

Create a wrapper for a K-value interpolator to be used with K-value flash.

The main purpose of this wrapper is to transform the general flash cond
NamedTuple into the right arguments for multi-linear interpolation.
"""
function KValueWrapper(K::T; dependence::Symbol = :pT) where T
    choices = (:p, :T, :pT, :pTz)
    dependence in choices || throw(ArgumentError("Bad input: dependence = $dependence was not in $choices"))
    return KValueWrapper{T, dependence}(K)
end

function (k::KValueWrapper{<:Any, :p})(cond::NamedTuple)
    val = k.K(cond.p)
    return val
end

function (k::KValueWrapper{<:Any, :T})(cond::NamedTuple)
    val = k.K(cond.T)
    return val
end

function (k::KValueWrapper{<:Any, :pT})(cond::NamedTuple)
    val = k.K(cond.p, cond.T)
    return val
end

function (k::KValueWrapper{<:Any, :T})(cond::NamedTuple)
    val = k.K(cond.p, cond.T, cond.z)
    return val
end
