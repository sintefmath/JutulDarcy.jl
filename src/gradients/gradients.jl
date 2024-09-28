include("objectives.jl")

"""
    result, sens = reservoir_sensitivities(case::JutulCase, objective::Function; sim_arg = NamedTuple(), kwarg...)

Simulate a case and calculate parameter sensitivities with respect to an
objective function on the form:

    obj(model, state, dt_n, n, forces_for_step_n)

The objective is summed up for all steps.

"""
function reservoir_sensitivities(case::JutulCase, objective::Function; sim_arg = NamedTuple(), kwarg...)
    result = simulate_reservoir(case; sim_arg...)
    (result, reservoir_sensitivities(case, result, objective; kwarg...))
end

"""
    reservoir_sensitivities(case::JutulCase, rsr::ReservoirSimResult, objective::Function; kwarg...)

Calculate parameter sensitivities with respect to an objective function on the
form for a case and a simulation result from that case. The objective function
is on the form:

    obj(model, state, dt_n, n, forces_for_step_n)

The objective is summed up for all steps.

    $(SIGNATURES)
"""
function reservoir_sensitivities(case::JutulCase, rsr::ReservoirSimResult, objective::Function; kwarg...)
    return reservoir_sensitivities(case, rsr.result, objective; kwarg...)
end

function reservoir_sensitivities(case::JutulCase, result::Jutul.SimResult, obj;
        include_parameters = false,
        include_state0 = false,
        adjoint_arg = (include_state0 = include_state0, ),
        kwarg...)
    sens = solve_adjoint_sensitivities(case, result, obj; adjoint_arg...)
    rmodel = reservoir_model(case.model)
    if haskey(sens, :Reservoir)
        sens = sens[:Reservoir]
    end
    data_domain_with_gradients = Jutul.data_domain_to_parameters_gradient(rmodel, sens; kwarg...)
    vars_to_add = []
    if include_parameters
        push!(vars_to_add, Jutul.get_parameters(rmodel))
    end
    if include_state0
        push!(vars_to_add, Jutul.get_primary_variables(rmodel))
    end
    for vars in vars_to_add
        for (k, pdef) in pairs(vars)
            u = associated_entity(pdef)
            data_domain_with_gradients[k, u] = sens[k]
        end
    end
    return data_domain_with_gradients
end
