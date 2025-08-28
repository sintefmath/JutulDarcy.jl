
import JutulDarcy: Tracers, get_phases,model_or_domain_is_well,phase_indices,upwind,face_average,gradient

using Jutul
import JutulDarcy: PhaseVariables


struct SurfactantTracer{N} <: Tracers.AbstractTracer
    phase_indices::NTuple{N, Int64}
    k_d::Float64   #Cg=k_d*Cl, where Cl is the surfactant concentration in liquid
    Concentration_limit::Union{String, Missing}
    C_max::Float64 #max tracer concentration, for gas or liquid
    D::Union{NamedTuple{(:gas, :liquid), Tuple{Float64, Float64}}, Nothing}  # Diffusivity
end


# C_max::Float64 
# This parameter is defined to prevent excessive accumulation of surfactant in the gas/liquid phase 
# For example: When injecting both gas and gas-soluable surfactant with a small k_d, it could lead to unreasonably high surfactant concentrations in liquid phase.
#
# Surfactant Tracer Transport Model
#
# This implementation models the transport of a surfactant tracer in a multiphase flow system.
# The surfactant distributes between the gas and liquid phases according to a partition coefficient k_d:
#
#     Cg = kd * Cl
#
# where:
# - Cg is the surfactant concentration in the gas phase.
# - Cl is the surfactant concentration in the liquid phase.
# - kd is the partition coefficient
#
# The effective tracer concentration in a given control volume is defined as:
#
#     Ceff = (Sl * ρl * Cl + Sg * ρg * Cg) / (Sl * ρl + Sg * ρg)
#
# The transport equation is given by:
# ∂(ϕ * Sl * ρl * Cl + ϕ * Sg * ρg * Cg) / ∂t + ∇ ⋅ (ul * ρl * Cl + ug * ρg * Cg) = 0
#
# Thus:
#  ∂(ϕ * Ceff *(Sl * ρl + Sg * ρg)) / ∂t + ∇ ⋅ (ul * ρl * Cl + ug * ρg * Cg) = 0
#
# - If the maximum surfactant tracer concentration is constrained in the liquid phase, it is computed as:
#    Cl = min((Ceff * (Sg * ρg + Sl * ρl)) / (k_d * Sg * ρg + Sl * ρl), C_max)
#    Cg = (Ceff * (Sg * ρg + Sl * ρl) - Sl * ρl * Cl) / (Sg * ρg)
#
# - If the maximum surfactant tracer concentration is constrained in the gas phase, it is computed as:
#    Cg = min((Ceff * (Sg * ρg + Sl * ρl) * k_d) / (k_d * Sg * ρg + Sl * ρl), C_max)
#    Cl = (Ceff * (Sg * ρg + Sl * ρl) - Sg * ρg * Cg) / (Sl * ρl)

function SurfactantTracer(system::MultiPhaseSystem, target_phases::AbstractArray;k_d=1.0,Concentration_limit=missing,C_max=Inf64,D=nothing)
    if !ismissing(Concentration_limit) && !(Concentration_limit in ("Liquid", "Gas"))
        error("Invalid Concentration_limit. Allowed values are Liquid or Gas")
    end
    if C_max < Inf && ismissing(Concentration_limit)
        error("C_max is finite ($C_max) but no Concentration_limit is specified. Use Liquid or Gas.")
    end
    nph = number_of_phases(system)
    N = length(target_phases)
    vals = Int[]
    for (p, target) in enumerate(target_phases)
        if target isa Int
            target <= nph || error("Phase index $target at position $p out of bounds")
            push!(vals, target)
        else
            ix = 0
            for (i, ph) in enumerate(get_phases(system))
                if ph == target
                    ix = i
                    break
                end
            end
            ix > 0 || error("Phase $target at position $p not found in system")
            push!(vals, ix)
        end
    end
    return SurfactantTracer{N}(tuple(vals...),k_d,Concentration_limit,C_max,D)
end

function SurfactantTracer(system::MultiPhaseSystem;k_d=1.0,Concentration_limit=missing,C_max=Inf64,D=nothing)
    if !ismissing(Concentration_limit) && !(Concentration_limit in ("Liquid", "Gas"))
        error("Invalid Concentration_limit. Allowed values are String Liquid or Gas.")
    end
    if C_max < Inf && ismissing(Concentration_limit)
        error("C_max is finite ($C_max) but no Concentration_limit is specified. Use Liquid or Gas.")
    end
    N = number_of_phases(system)
    return SurfactantTracer{N}(tuple((1:N)...),k_d,Concentration_limit,C_max,D)
end

Tracers.tracer_phase_indices(t::SurfactantTracer) = t.phase_indices


struct LiquidSurfConcentration<: Jutul.ScalarVariable
    tracer_ix::Int
end

Jutul.minimum_value(t::LiquidSurfConcentration) = 0.0

function Jutul.update_secondary_variable!(cl, def::LiquidSurfConcentration, model, state, ix)
    TracerConcentrations = state.TracerConcentrations
    surf_tracer=model.equations[:tracers].flux_type.tracers[def.tracer_ix]
    C_max= surf_tracer.C_max # maximum tracer concentration in liquid/gas phase
    gas_index = findfirst(ph -> ph == VaporPhase(), get_phases(model.system))
    k_d=surf_tracer.k_d
    Concentration_limit=surf_tracer.Concentration_limit
    cg=state.GasSurfConcentration
    for cell in ix
        ceff = TracerConcentrations[def.tracer_ix, cell]
        (liquid_mass, gas_mass)=compute_phase_masses(surf_tracer, state, cell, gas_index)
        if ismissing(Concentration_limit) || Concentration_limit == "Liquid"
            cl[cell]=min(ceff*(gas_mass+liquid_mass)/(k_d*gas_mass+liquid_mass),C_max)
        else
            if liquid_mass > 0
                cl[cell]=(ceff*(gas_mass+liquid_mass)-cg[cell]*gas_mass)/liquid_mass
            else
                cl[cell]=0.0
            end
        end
    end
    return cl
end


struct GasSurfConcentration <: Jutul.ScalarVariable
    tracer_ix::Int
end

Jutul.minimum_value(t::GasSurfConcentration) = 0.0


function Jutul.update_secondary_variable!(cg, def::GasSurfConcentration, model, state, ix)
    TracerConcentrations = state.TracerConcentrations

    surf_tracer=model.equations[:tracers].flux_type.tracers[def.tracer_ix]
    C_max= surf_tracer.C_max # Wanted max tracerconcentration in liquid phase
    gas_index = findfirst(ph -> ph == VaporPhase(), get_phases(model.system))
    k_d=surf_tracer.k_d
    cl=state.LiquidSurfConcentration
    Concentration_limit=surf_tracer.Concentration_limit
    for cell in ix
        ceff = TracerConcentrations[def.tracer_ix, cell]
        (liquid_mass, gas_mass)=compute_phase_masses(surf_tracer, state, cell, gas_index)
        if !ismissing(Concentration_limit) && Concentration_limit == "Gas"
            cg[cell]=min((ceff*(gas_mass+liquid_mass)*k_d)/(k_d*gas_mass+liquid_mass),C_max)
        else
            if gas_mass > 0
                cg[cell] = (ceff*(gas_mass+liquid_mass) - liquid_mass*cl[cell]) / gas_mass
            else
                cg[cell] = 0.0  
            end
        end
    end
    return cg
end


function compute_phase_masses(surf_tracer, state, cell, gas_index)
    liquid_mass = 0.0
    gas_mass = 0.0
    for phase in Tracers.tracer_phase_indices(surf_tracer)
        S = state.Saturations[phase, cell]
        den = state.PhaseMassDensities[phase, cell]
        if phase == gas_index
            gas_mass = S * den
        else
            liquid_mass += S * den
        end
    end
    return (liquid_mass, gas_mass)
end





function tracer_total_mass_flux(tracer::SurfactantTracer, model, state, phase_mass_fluxes, index, upw,kgrad,face)
    D = tracer.D
    cg=state.GasSurfConcentration
    cl=state.LiquidSurfConcentration
    gas_index = findfirst(ph -> ph == VaporPhase(), get_phases(model.system))
    v = zero(eltype(cg))
    for phase in tracer_phase_indices(tracer)
        q_ph = phase_mass_fluxes[phase]
        if phase == gas_index
            F= cell -> cg[cell]
                else
            F= cell -> cl[cell]
        end
        C_iface = upwind(upw, F, q_ph)
        v += C_iface * q_ph
        v = add_diffusive_tracer_flux(v,state,D, face, kgrad, phase,cg,cl,gas_index)
    end
    return v
end

#diffusion

function add_diffusive_tracer_flux(v,state,D::Nothing, face, kgrad, phase,cg,cl,gas_index)
    return v
end


function add_diffusive_tracer_flux(v,state,D, face, kgrad, phase,cg,cl,gas_index)
    
    if phase == gas_index
        C = cg
        D_phase = D.gas
    else
        C = cl
        D_phase = D.liquid
    end

    diff_mass = phase_diffused_mass(D_phase, state, face, kgrad,phase)
    dC = gradient(C, kgrad)
    v += diff_mass*dC
    return v
end

function phase_diffused_mass(D_phase, state, face, kgrad,phase)
    S = state.Saturations
    ρ=state.PhaseMassDensities
    den = cell -> @inbounds ρ[phase, cell]
    left, right = Jutul.cell_pair(kgrad)
    Sa = min(S[phase, left], S[phase, right])
    return -D_phase*Sa*face_average(den, kgrad)
end

    

#CMG-STARS surfactant model
#Fmfactor should works on gas relative permeability
#But in this model Fmfactor works on gas viscosity


struct FmFactor{V}<:Jutul.ScalarVariable
    epdry::V
    fmdry::V
    fmsurf::V   # kg/m3 or g/l, critical surfactant concentration in liquid
    epsurf::V
    fommb::V
    #swr::V #residual water saturation
end


function FmFactor(; epdry=500.0, fmdry=0.25, fmsurf=2.5, epsurf=1.0, fommb=500.0)
    return FmFactor(epdry, fmdry, fmsurf, epsurf, fommb)
end


function Jutul.default_value(model, v::FmFactor)
    return 1.0
end 


@jutul_secondary function update_Fmfactor!(fm, Fm_def::FmFactor, model, ClVolumetricConcentration, Saturations, cells)
    epdry, fmdry, fmsurf, epsurf, fommb = Fm_def.epdry, Fm_def.fmdry, Fm_def.fmsurf, Fm_def.epsurf, Fm_def.fommb
    phases = get_phases(model.system)
    gas_index = findfirst(ph -> ph == VaporPhase(), phases)

    if length(phases) != 2 || gas_index === nothing
        error("This Fmfactor only work for two-phase system with VaporPhase present")
    end
    
    liquid_index = 3 - gas_index
    
    for i in cells
        Csw = ClVolumetricConcentration[i]
        sl = Saturations[liquid_index, i]
        F_dry = get_F_dry(sl, epdry, fmdry)
        F_surf = get_F_surfactant(Csw, fmsurf, epsurf)
        fm[i] = 1 / (1 + fommb * F_dry * F_surf)
    end
    
    return fm
end



function get_F_dry(sl::T, epdry::Real, fmdry::Real) where T
    sl = clamp(sl, zero(T), one(T))
    return 0.5+atan((epdry*(sl-fmdry)))/pi
end 



function get_F_surfactant(C_sw::T, fmsurf::Real, epsurf::Real) where T
    if C_sw<=fmsurf
        F_surf=(C_sw/fmsurf)^epsurf
    else
        F_surf=1.0
    end
    return F_surf
end





struct ClVolumetricConcentration<: Jutul.ScalarVariable end


Jutul.minimum_value(t::ClVolumetricConcentration) = 0.0

function Jutul.default_value(model, v::ClVolumetricConcentration)
    return 1.0
end 

function Jutul.update_secondary_variable!(cl_vol, def::ClVolumetricConcentration, model,state,ix)
    cl=state.LiquidSurfConcentration
    dens=state.PhaseMassDensities
    gas_index = findfirst(ph -> ph == VaporPhase(), get_phases(model.system))
    liquid_index=3-gas_index
    for cell in ix
        cl_vol[cell]=cl[cell]*dens[liquid_index,cell]
        #This function converts the mass fraction (kg/kg) of surfactant in the liquid phase to a volumetric concentration (g/L).
    end
    return cl_vol
end




####
# In CMG-STARS, `FmFactor` is used to modify the gas relative permeability:
#
#     Krg = Krg * FmFactor
#
# In this model, `FmFactor` is instead used to modify the gas viscosity:
#
#     gas_viscosity = gas_viscosity / FmFactor
#
# This adjustment affects the mobility of the gas phase in the two-phase Darcy flow equation:
#
#     v_α = - (k_rα / μ_α) * K ⋅ ∇p
#
# where:
# - v_α is the velocity of phase α,
# - k_rα is the relative permeability of phase α,
# - μ_α is the viscosity of phase α,
# - K is the absolute permeability tensor,
# - ∇p is the pressure gradient.
#
# By modifying the gas viscosity instead of its relative permeability, `FmFactor` still controls gas mobility but in a different manner.


struct SurfactantAdjustedViscosities <: JutulDarcy.PhaseVariables
end

Jutul.@jutul_secondary function update_gas_viscosity!(vals, def::SurfactantAdjustedViscosities, model, BasePhaseViscosities, FmFactor, ix)

    phases = get_phases(model.system)
    gas_index = findfirst(ph -> ph == VaporPhase(), phases)
    nph = size(vals, 1)

    for cell in ix
        for ph in 1:nph
            val_phase = BasePhaseViscosities[ph, cell]
            if ph == gas_index 
                val_phase /=FmFactor[cell]
            end
            vals[ph, cell] = val_phase
        end
    end
    return vals
end


function set_surfactant_model!(outer_model::MultiModel; is_well = false)
    rmodel = reservoir_model(outer_model)
    set_surfactant_model!(rmodel, is_well = false)
    for (mname, model) in pairs(outer_model.models)
        if JutulDarcy.model_or_domain_is_well(model)
            set_surfactant_model!(model, is_well = true)
        end
    end
    return outer_model
end


function set_surfactant_model!(model::SimulationModel;is_well=false)
    tracer_ix = findfirst(x -> isa(x, SurfactantTracer), model.equations[:tracers].flux_type.tracers)
    set_secondary_variables!(model,
    GasSurfConcentration = GasSurfConcentration(tracer_ix),
    LiquidSurfConcentration = LiquidSurfConcentration(tracer_ix),
    ClVolumetricConcentration=ClVolumetricConcentration()
    )

    if !is_well
        prms = Jutul.get_variables_by_type(model, :parameters)
        svars = Jutul.get_variables_by_type(model, :secondary)
        set_secondary_variables!(model,
        FmFactor=FmFactor(),)
        if haskey(prms, :PhaseViscosities)
            mu = prms[:PhaseViscosities]
            prms[:BasePhaseViscosities] = mu
            delete!(prms, :PhaseViscosities)
        else
            @assert haskey(svars, :PhaseViscosities) "PhaseViscosities not found in secondary variables or parameters, cannot setup foam-surfactant model."
            mu = svars[:PhaseViscosities]
            svars[:BasePhaseViscosities] = mu
            delete!(svars, :PhaseViscosities)
        end
        svars[:PhaseViscosities] = SurfactantAdjustedViscosities()
    end
    push!(model.output_variables, :LiquidSurfConcentration)
    push!(model.output_variables, :GasSurfConcentration)
    return model
end


function tracer_scale(model, tracer::SurfactantTracer)

    return 100.0
end


