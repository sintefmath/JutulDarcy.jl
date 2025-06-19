struct PressureReservoirFromWellFlowCT{T} <: Jutul.AdditiveCrossTerm
    parent::T
end

function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t,
        state_s, state0_s, 
        model_t, model_s,
        ct::PressureReservoirFromWellFlowCT,
        eq::PressureEquation, dt, ldisc = local_discretization(ct, i)
    )

    sys = model_t.system
    rhoS = reference_densities(sys)
    conn = cross_term_perforation_get_conn(ct.parent, i, state_s, state_t)

    ncomp = number_of_components(sys)
    out_full = zeros(eltype(out), ncomp)
    @assert !haskey(state_s, :MassFractions)
    q = multisegment_well_perforation_flux!(out_full, sys, state_t, state_s, rhoS, conn)
    out[1] = sum(out_full)
end


function Jutul.cross_term_entities(ct::PressureReservoirFromWellFlowCT, eq::PressureEquation, model)
    return ct.parent.reservoir_cells
end

function Jutul.cross_term_entities_source(ct::PressureReservoirFromWellFlowCT, eq::PressureEquation, model)
    return ct.parent.well_cells
end

Jutul.symmetry(::PressureReservoirFromWellFlowCT) = Jutul.CTSkewSymmetry()

# function Jutul.setup_cross_term_storage(
#         ct::PressureReservoirFromWellFlowCT,
#         eq_t::PressureEquation,
#         eq_s,
#         model_t,
#         model_s,
#         storage_t,
#         storage_s
#     )
#     eq_t_cons = eq_t.conservation

#     s = Jutul.setup_cross_term_storage(
#         ct.parent,
#         eq_t_cons,
#         eq_s,
#         model_t,
#         model_s,
#         storage_t,
#         storage_s
#         )
#     return s
# end


# function Jutul.align_to_jacobian!(
#     eq_s,
#     ct_p::PressureReservoirFromWellFlowCT,
#     jac,
#     model,
#     entity,
#     impact_t;
#     context = model.context,
#     positions = nothing,
#     equation_offset = 0,
#     variable_offset = 0,
#     diagonal = true,
#     number_of_entities_target = nothing,
#     kwarg...
#     )
#     ct = ct_p.parent
#     k = Jutul.entity_as_symbol(entity)
#     has_pos = !isnothing(positions)
#     if haskey(eq_s, k)
#         cache = eq_s[k]
#         if has_pos
#             # Align against other positions that is provided
#             pos = positions[k]
#         else
#             # Use positions from cache
#             pos = cache.jacobian_positions
#         end
#         # J = cache.variables
#         if isnothing(number_of_entities_target)
#             nt = cache.number_of_entities_target
#         else
#             nt = number_of_entities_target
#         end
#         I, J = Jutul.generic_cache_declare_pattern(cache, impact_t)

#         Nt, Ne, Np = Jutul.ad_dims(cache)
#         # Diagonal Reservoir wrt Reservoir: Ne = Np = 1
#         # Off-diagonal reservoir wrt. wells Ne = Np = 1
#         is_reservoir = model isa PressureModel
#         is_pressure_part = (is_reservoir && diagonal) || (!is_reservoir && !diagonal)
#         is_pressure_part = is_reservoir
#         if is_pressure_part
#             # Ne = Np = 1 # Assume 1 eq, 1 pressure atm
#         end
#         reservoir_diagonal = is_reservoir && diagonal
#         reservoir_off_diagonal = !is_reservoir && !diagonal
#         well_diagonal = !is_reservoir && diagonal
#         well_off_diagonal = is_reservoir && !diagonal

#         if reservoir_diagonal
#             Np = 1
#             Ne = 1
#         elseif reservoir_off_diagonal
#             Ne = 1
#         elseif well_off_diagonal
#             Np = 1
#         else
#             @assert well_diagonal
#         end
#         dims = (Nt, Ne, Np)


#         Jutul.injective_alignment!(
#             cache,
#             ct,
#             jac,
#             entity,
#             context;
#             pos = pos,
#             target_index = I,
#             source_index = J,
#             number_of_entities_source = cache.number_of_entities_source,
#             number_of_entities_target = nt,
#             target_offset = equation_offset,
#             source_offset = variable_offset,
#             dims = dims,
#             kwarg...
#         )
#     else
#         @warn "Did not find $k in $(keys(eq_s))"
#     end
# end
