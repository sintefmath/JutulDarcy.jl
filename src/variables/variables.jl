include("capillary.jl")
include("relperm.jl")
include("pvt.jl")
include("utils.jl")
include("viscosity.jl")

degrees_of_freedom_per_entity(model, sf::PhaseVariables) = number_of_phases(model.system)

# Generic version
degrees_of_freedom_per_entity(model, sf::ComponentVariables) = number_of_components(model.system)

# Single-phase specialization
degrees_of_freedom_per_entity(model::SimulationModel{D, S}, sf::ComponentVariables) where {D, S<:SinglePhaseSystem} = 1

# Immiscible specialization
degrees_of_freedom_per_entity(model::SimulationModel{D, S}, sf::ComponentVariables) where {D, S<:ImmiscibleSystem} = number_of_phases(model.system)

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
    is_well = physical_representation(domain) isa WellDomain
    is_multiphase = !isa(system, SinglePhaseSystem)
    is_bo = system isa BlackOilSystem
    if is_multiphase && !is_well
        S[:RelativePermeabilities] = BrooksCoreyRelativePermeabilities(system)
    end
    if !is_well
        S[:PhaseMobilities] = PhaseMobilities()
    end
    if !is_bo && !is_well
        S[:PhaseMassMobilities] = PhaseMassMobilities()
    end
end

function select_default_darcy_parameters!(prm, domain, system::SinglePhaseSystem, formulation)
    prm[:PhaseViscosities] = PhaseViscosities()
    prm[:FluidVolume] = FluidVolume()
    prm[:RelativePermeabilities] = BrooksCoreyRelativePermeabilities(system)
    prm[:Saturations] = Saturations()
end

function select_default_darcy_parameters!(prm, domain, system::ImmiscibleSystem, formulation)
    add_connate_water_if_aqueous_present!(prm, domain, system)
    prm[:PhaseViscosities] = PhaseViscosities()
    prm[:FluidVolume] = FluidVolume()
    if number_of_phases(system) == 1
        prm[:Saturations] = Saturations()
    end
end

function select_default_darcy_parameters!(prm, domain, system::MultiPhaseSystem, formulation)
    add_connate_water_if_aqueous_present!(prm, domain, system)
    prm[:FluidVolume] = FluidVolume()
    if number_of_phases(system) == 1
        prm[:Saturations] = Saturations()
    end
end

function add_connate_water_if_aqueous_present!(prm, domain, system)
    is_well = model_or_domain_is_well(domain)
    if !is_well
        if AqueousPhase() in get_phases(system)
            prm[:ConnateWater] = ConnateWater()
        end
    end
end

function select_minimum_output_variables!(out, system::MultiPhaseSystem, model)
    push!(out, :TotalMasses)
end

function model_or_domain_is_well(m)
    if m isa SimulationModel
        domain = m.domain
    else
        domain = m
    end
    if domain isa DiscretizedDomain
        is_well = physical_representation(domain) isa WellDomain
    else
        is_well = false
    end
    return is_well
end
