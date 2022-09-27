export PhaseMassDensities, ConstantCompressibilityDensities
export BrooksCoreyRelPerm, TabulatedRelPermSimple

include("capillary.jl")
include("relperm.jl")
include("pvt.jl")
include("utils.jl")
include("viscosity.jl")

degrees_of_freedom_per_entity(model, sf::PhaseVariables) = number_of_phases(model.system)

# Single-phase specialization
degrees_of_freedom_per_entity(model::SimulationModel{D, S}, sf::ComponentVariable) where {D, S<:SinglePhaseSystem} = 1

# Immiscible specialization
degrees_of_freedom_per_entity(model::SimulationModel{D, S}, sf::ComponentVariable) where {D, S<:ImmiscibleSystem} = number_of_phases(model.system)

function select_secondary_variables!(S, system::MultiPhaseSystem, model)
    select_default_darcy_secondary_variables!(S, model.domain, system, model.formulation)
end

function select_parameters!(S, system::MultiPhaseSystem, model)
    select_default_darcy_parameters!(S, model.domain, system, model.formulation)
end

function select_default_darcy_secondary_variables!(S, domain, system, formulation)
    nph = number_of_phases(system)
    S[:PhaseMassDensities] = ConstantCompressibilityDensities(nph)
    S[:TotalMasses] = TotalMasses()
    if !(isa(system, SinglePhaseSystem) || isa(domain.grid, WellGrid))
        S[:RelativePermeabilities] = BrooksCoreyRelPerm(system)
    end
end

function select_default_darcy_parameters!(prm, domain, system::SinglePhaseSystem, formulation)
    prm[:PhaseViscosities] = PhaseViscosities()
    prm[:FluidVolume] = FluidVolume()
    prm[:RelativePermeabilities] = BrooksCoreyRelPerm(system)
    prm[:Saturations] = Saturations()
end

function select_default_darcy_parameters!(prm, domain, system::ImmiscibleSystem, formulation)
    prm[:PhaseViscosities] = PhaseViscosities()
    prm[:FluidVolume] = FluidVolume()
end

function select_default_darcy_parameters!(prm, domain, system::MultiPhaseSystem, formulation)
    prm[:FluidVolume] = FluidVolume()
end

function select_minimum_output_variables!(out, system::MultiPhaseSystem, model)
    push!(out, :TotalMasses)
end
