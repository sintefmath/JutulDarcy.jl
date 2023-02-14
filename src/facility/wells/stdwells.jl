const SimpleWellDomain = DiscretizedDomain{<:SimpleWell}
const SimpleWellFlowModel = SimulationModel{<:SimpleWellDomain, <:MultiPhaseSystem}

struct WellMassFractions <: FractionVariables end
Jutul.need_default_primary(model, ::WellMassFractions) = false
function Jutul.default_values(model, mf::WellMassFractions)
    nc = values_per_entity(model, mf)
    return [1/nc for _ in 1:nc]
end

struct SimpleWellEquation <: JutulEquation end

Jutul.local_discretization(::SimpleWellEquation, i) = nothing

function Jutul.number_of_equations_per_entity(model::SimpleWellFlowModel, ::SimpleWellEquation)
    sys = model.system
    return number_of_components(sys)
end

function values_per_entity(model, v::WellMassFractions)
    sys = model.system
    return number_of_components(sys)
end

function select_primary_variables!(model::SimpleWellFlowModel)
    pvars = model.primary_variables
    pvars[:Pressure] = Pressure()
    pvars[:MassFractions] = WellMassFractions()
end

function select_secondary_variables!(model::SimpleWellFlowModel)

end

function select_equations!(model::SimpleWellFlowModel)
    model.equations[:mass_conservation] = SimpleWellEquation()
end

function select_parameters!(model::SimpleWellFlowModel)
    prm = model.parameters
    prm[:WellIndices] = WellIndices()
    prm[:PerforationGravityDifference] = PerforationGravityDifference()
end

function select_minimum_output_variables!(model::SimpleWellFlowModel)
    outputs = model.output_variables
    for k in keys(model.primary_variables)
        push!(outputs, k)
    end
    return model
end

