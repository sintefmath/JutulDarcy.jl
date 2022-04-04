
get_components(sys::MultiComponentSystem) = sys.components
number_of_components(sys::MultiComponentSystem) = length(get_components(sys))

liquid_phase_index(sys::MultiPhaseCompositionalSystemLV) = phase_index(sys, LiquidPhase())
vapor_phase_index(sys::MultiPhaseCompositionalSystemLV) = phase_index(sys, VaporPhase())
other_phase_index(sys::MultiPhaseCompositionalSystemLV{E, T, O}) where {E, T, O} = phase_index(sys, O())

phase_index(sys, phase) = only(findfirst(isequal(phase), sys.phases))
has_other_phase(sys) = true
has_other_phase(sys::MultiPhaseCompositionalSystemLV{E, T, O}) where {E, T, O<:Nothing} = false

phase_indices(sys::MultiPhaseCompositionalSystemLV{E, T, O}) where {E, T, O<:Nothing} = (liquid_phase_index(sys), vapor_phase_index(sys))
phase_indices(sys) = (other_phase_index(sys), liquid_phase_index(sys), vapor_phase_index(sys))
