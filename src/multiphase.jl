
# Abstract multiphase system

get_phases(sys::MultiPhaseSystem) = sys.phases
get_phases(s::Symbol) = get_phases(Val(s))

@inline number_of_phases(sys::MultiPhaseSystem) = length(get_phases(sys))
@inline reference_densities(sys::MultiPhaseSystem) = sys.rho_ref
@inline reference_densities(sys::CompositeSystem) = reference_densities(flow_system(sys))

flow_system(sys::MultiPhaseSystem) = sys
flow_system(sys::CompositeSystem) = sys.systems.flow

number_of_components(sys::ImmiscibleSystem) = number_of_phases(sys)
number_of_components(sys::CompositeSystem) = number_of_components(flow_system(sys))

"""
    component_names(sys)

Get a list of the component names (as Strings)
"""
component_names(sys::Union{SinglePhaseSystem, ImmiscibleSystem}) = phase_names(sys)
component_names(sys::CompositeSystem) = phase_names(sys.systems.flow)

phase_names(system) = phase_name.(get_phases(system))

get_phases(sys::SinglePhaseSystem) = (sys.phase, )

phase_indices(sys::SinglePhaseSystem) = 1
phase_indices(sys::ImmiscibleSystem) = tuple(eachindex(sys.phases)...)

number_of_phases(::SinglePhaseSystem) = 1
number_of_phases(sys::CompositeSystem) = number_of_phases(sys.systems.flow)

"""
    eachphase(sys::MultiPhaseSystem)

Get tuple of the phases (1,...,Nph). Convenient when you want statically
known compile time iteration, for example by use of `map`.
"""
eachphase(sys::MultiPhaseSystem) = tuple(eachindex(sys.phases)...)
eachphase(sys::SinglePhaseSystem) = (1,)

"""
    AqueousPhase()

`AbstractPhase` subtype for water-like phases.
"""
struct AqueousPhase <: AbstractPhase end
phase_name(::AqueousPhase) = "Aqueous"

"""
    LiquidPhase()

`AbstractPhase` subtype for liquid-like phases.
"""
struct LiquidPhase <: AbstractPhase end
phase_name(::LiquidPhase) = "Liquid"

"""
    VaporPhase()

`AbstractPhase` subtype for vapor or gaseous phases.
"""
struct VaporPhase <: AbstractPhase end
phase_name(::VaporPhase) = "Vapor"

## Main implementation
# Primary variable logic

const DEFAULT_MINIMUM_PRESSURE = 10000.0 # 0.1 MPa

struct Pressure <: ScalarVariable
    max_abs::Union{Float64, Nothing}
    max_rel::Union{Float64, Nothing}
    minimum_pressure::Float64
    maximum_pressure::Float64
    scale::Float64
    function Pressure(max_abs, max_rel, minimum_pressure, maximum_pressure, scale)
        if !isnothing(max_abs)
            @assert max_abs > 0.0 "Maximum absolute pressure change was $max_abs, must be positive"
        end
        if !isnothing(max_rel)
            @assert max_rel > 0.0 "Maximum relative pressure change was $max_rel, must be positive"
        end
        @assert minimum_pressure < maximum_pressure "Maximum pressure $maximum_pressure must be larger than minimum pressure $minimum_pressure"
        @assert scale > 0.0 "Pressure scale must be positive Float64, was $scale"
        new(max_abs, max_rel, minimum_pressure, maximum_pressure, scale)
    end
end

"""
    Pressure(; max_abs = nothing, max_rel = 0.2, scale = si_unit(:bar), maximum = Inf, minimum = DEFAULT_MINIMUM_PRESSURE)

Pressure variable definition. `max_abs`/`max_rel` maximum allowable
absolute/relative change over a Newton iteration, `scale` is a "typical" value
used to regularize the linear system, `maximum` the largest possible value and
`minimum` the smallest.
"""
function Pressure(; max_abs = nothing, max_rel = 0.2, scale = si_unit(:bar), maximum = Inf, minimum = DEFAULT_MINIMUM_PRESSURE)
    return Pressure(max_abs, max_rel, minimum, maximum, scale)
end

Jutul.variable_scale(p::Pressure) = p.scale
absolute_increment_limit(p::Pressure) = p.max_abs
relative_increment_limit(p::Pressure) = p.max_rel
maximum_value(p::Pressure) = p.maximum_pressure
minimum_value(p::Pressure) = p.minimum_pressure

struct Saturations <: FractionVariables
    ds_max::Float64
end

"""
    Saturations(; ds_max = 0.2)

Saturations as primary variable. `ds_max` controls maximum allowable saturation
change between two successive Newton iterations.
"""
function Saturations(;ds_max = 0.2)
    Saturations(ds_max)
end

function default_value(model, ::Saturations)
    nph = number_of_phases(model.system)
    return 1.0/nph
end

values_per_entity(model, v::Saturations) = number_of_phases(model.system)
absolute_increment_limit(s::Saturations) = s.ds_max

function initialize_primary_variable_ad!(state, model, pvar::Saturations, state_symbol, npartials; kwarg...)
    nph = values_per_entity(model, pvar)
    v = state[state_symbol]
    state[state_symbol] = Jutul.unit_sum_init(v, model, npartials, nph; kwarg...)
    return state
end

Jutul.parameter_is_differentiable(::Saturations, model) = false

"""
    ConnateWater()

Parameter for connate water per cell. Used in some three-phase relative
permeability evaluators.
"""
struct ConnateWater <: ScalarVariable end

function Jutul.default_values(model, ::ConnateWater)
    nc = number_of_cells(model.domain)
    relperm = Jutul.get_variable(model, :RelativePermeabilities)
    swcon = zeros(nc)
    if hasproperty(relperm, :regions)
        kr = relperm[:w]
        for i in 1:nc
            reg = JutulDarcy.region(relperm.regions, i)
            kr_i = JutulDarcy.table_by_region(kr, reg)
            swcon[i] = kr_i.connate
        end
    end
    return swcon
end

"""
    TotalMasses()

Variable that defines total component masses in each cell of the domain.
"""
struct TotalMasses <: VectorVariables end

function degrees_of_freedom_per_entity(model::SimulationModel{G, S}, v::TotalMasses) where {G<:Any, S<:MultiPhaseSystem}
    number_of_phases(model.system)
end

struct PhaseMassMobilities <: PhaseVariables end
struct PhaseMobilities <: PhaseVariables end

@inline function minimum_value(::TotalMasses)
    0.0
end

"""
    TotalMasses()

Variable that defines total mass of all components in each cell of the domain.
"""
struct TotalMass <: ScalarVariable end
@inline function minimum_value(::TotalMass)
    return 0.0
end

"""
    Transmissibilities()

Variable/parameter used to define the cell-to-cell transmissibilities when using
a two-point flux approximation scheme.
"""
struct Transmissibilities <: ScalarVariable end
Jutul.variable_scale(::Transmissibilities) = 1e-10
Jutul.minimum_value(::Transmissibilities) = 0.0

Jutul.associated_entity(::Transmissibilities) = Faces()
function Jutul.default_values(model, ::Transmissibilities)
    return physical_representation(model.domain).trans
end

function Jutul.default_parameter_values(data_domain, model, param::Transmissibilities, symb)
    if haskey(data_domain, :transmissibilities, Faces())
        # This takes precedence
        T = copy(data_domain[:transmissibilities])
    elseif haskey(data_domain, :permeability, Cells())
        T = reservoir_transmissibility(data_domain)
    else
        error(":permeability or :transmissibilities symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    if eltype(T)<:AbstractFloat
        replace_bad_trans!(T, symb)
    end
    return T
end

function replace_bad_trans!(T, symb; replace = 0.0)
    for (F, descr) in [(x -> x < 0, "negative"), (x -> !isfinite(x), "non-finite")]
        if any(F, T)
            c = count(F, T)
            @warn "Parameter initialization for $symb: $c $descr values detected out of $(length(T)) total. Replacing with $replace"
            for i in eachindex(T)
                if F(T[i])
                    T[i] = replace
                end
            end
        end
    end
    return T
end

struct Diffusivities <: VectorVariables end
Jutul.variable_scale(::Diffusivities) = 1e-10
Jutul.minimum_value(::Diffusivities) = 0.0

Jutul.associated_entity(::Diffusivities) = Faces()
Jutul.values_per_entity(model, ::Diffusivities) = number_of_phases(model.system)

function Jutul.default_parameter_values(data_domain, model, param::Diffusivities, symb)
    nf = number_of_faces(model.domain)
    nph = number_of_phases(model.system)
    if haskey(data_domain, :diffusivities, Faces())
        # This takes precedence
        T = data_domain[:diffusivities]
    elseif haskey(data_domain, :diffusion, Cells())
        T = zeros(nph, nf)
        ϕ = data_domain[:porosity]
        D = data_domain[:diffusion]
        U = ϕ'.*D
        g = physical_representation(data_domain)
        if U isa AbstractVector
            T_i = compute_face_trans(g, U)
            for i in 1:nph
                T[i, :] .= T_i
            end
        else
            for i in 1:nph
                T_i = compute_face_trans(g, U[i, :])
                T[i, :] .= T_i
            end
        end
        if any(x -> x < 0, T)
            c = count(x -> x < 0, T)
            @warn "$c negative diffusivities detected."
        end
    else
        error(":diffusion or :diffusivities symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    @assert size(T) == (nph, nf)
    return T
end

"""
    TwoPointGravityDifference()

Parameter representing the difference in gravity on an instance of `Faces`
between two `Cells`. If the phase flux is written as

``v = - K \\nabla (p + \\rho g \\nabla z)``

this term represent the discretized analogue of ``\\rho g \\nabla z``.
"""
struct TwoPointGravityDifference <: ScalarVariable end

Jutul.associated_entity(::TwoPointGravityDifference) = Faces()
function Jutul.default_values(model, ::TwoPointGravityDifference)
    return physical_representation(model.domain).gdz
end

function Jutul.default_parameter_values(data_domain, model, param::TwoPointGravityDifference, symb)
    has_gdz = haskey(data_domain, :gdz, Faces())
    has_cell_centroids = haskey(data_domain, :cell_centroids, Cells())
    has_face_neighbors = haskey(data_domain, :neighbors, Faces())
    if has_gdz
        # This takes precedence
        gdz = data_domain[:gdz]
    elseif has_cell_centroids && has_face_neighbors
        N = data_domain[:neighbors]
        cc = data_domain[:cell_centroids]
        if size(cc, 1) == 3
            z = vec(cc[3, :])
            gdz = compute_face_gdz(N, z)
        else
            nf = size(N, 2)
            gdz = zeros(nf)
        end
    else
        error(":gdz or :neighbors + :cell_centroids symbols must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return gdz
end

# Selection of variables
function select_primary_variables!(S, ::SinglePhaseSystem, model::SimulationModel)
    S[:Pressure] = Pressure()
end

function select_primary_variables!(S, sys::ImmiscibleSystem, model::SimulationModel)
    S[:Pressure] = Pressure()
    if number_of_phases(sys) > 1
        S[:Saturations] = Saturations()
    end
end

function select_equations!(eqs, sys::MultiPhaseSystem, model::SimulationModel)
    fdisc = model.domain.discretizations.mass_flow
    nc = number_of_components(sys)
    eqs[:mass_conservation] = ConservationLaw(fdisc, :TotalMasses, nc)
end

function select_parameters!(prm, disc::D, model::DarcyFlowModel) where D<:Union{TwoPointPotentialFlowHardCoded, Jutul.PotentialFlow}
    prm[:Transmissibilities] = Transmissibilities()
    prm[:TwoPointGravityDifference] = TwoPointGravityDifference()
end

number_of_equations_per_entity(system::MultiPhaseSystem, e::ConservationLaw) = number_of_components(system)
number_of_equations_per_entity(system::SinglePhaseSystem, e::ConservationLaw) = 1

function pore_volume(data_domain::DataDomain; throw = true)
    if haskey(data_domain, :pore_volume, Cells())
        pv = data_domain[:pore_volume]
    elseif haskey(data_domain, :volumes, Cells())
        vol = copy(data_domain[:volumes])
        ntg = poro = pvmult = 1.0
        if haskey(data_domain, :porosity, Cells())
            poro = data_domain[:porosity]
        end
        if haskey(data_domain, :net_to_gross, Cells())
            ntg = data_domain[:net_to_gross]
        end
        if haskey(data_domain, :pore_volume_multiplier, Cells())
            pvmult = data_domain[:pore_volume_multiplier]
        end
        pv = @. vol * poro * ntg * pvmult
        if haskey(data_domain, :pore_volume_override, Cells())
            for (i, v) in enumerate(data_domain[:pore_volume_override])
                if isfinite(v)
                    pv[i] = v
                end
            end
        end
    else
        pv = missing
        if throw
            error("Neither pair :volumes and :porosity or :pore_volume found in domain.")
        end
    end

    return pv
end

pore_volume(model::MultiModel, parameters) = pore_volume(reservoir_model(model), parameters[:Reservoir])
pore_volume(model::SimulationModel, parameters) = fluid_volume(model, parameters)
fluid_volume(model, parameters) = parameters[:FluidVolume]
domain_fluid_volume(g) = missing

function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel{D, S}, eq::ConservationLaw, eq_s, force::V, time) where {V <: AbstractVector{SourceTerm{I, F, T}}, D, S<:MultiPhaseSystem} where {I, F, T}
    state = storage.state
    if haskey(state, :RelativePermeabilities)
        kr = state.RelativePermeabilities
    else
        kr = 1.0
    end
    mu = state.PhaseViscosities
    rhoS = reference_densities(model.system)
    insert_phase_sources!(acc, model, kr, mu, rhoS, force)
end

function local_mobility(kr::Real, mu, ph, c)
    return kr./mu[ph, c]
end

function local_mobility(kr, mu, ph, c)
    return kr[ph, c]./mu[ph, c]
end

function phase_source(c, src, rhoS, kr, mu, ph)
    v = src.value
    # c = src.cell
    if v > 0
        q = in_phase_source(src, v, c, kr, mu, ph)
    else
        q = out_phase_source(src, v, c, kr, mu, ph)
    end
    return rhoS*q
end


function in_phase_source(src, v, c, kr, mu, ph)
    f = src.fractional_flow[ph]
    return v*f
end

function out_phase_source(src, v, c, kr, mu, ph)
    mobT = 0
    mob = 0
    for i = 1:size(mu, 1)
        mi = local_mobility(kr, mu, i, c)
        mobT += mi
        if ph == i
            mob = mi
        end
    end
    f = mob/mobT
    return v*f
end

function insert_phase_sources!(acc, model, kr, mu, rhoS, sources)
    nph = size(acc, 1)
    M = global_map(model.domain)
    for src in sources
        c = Jutul.full_cell(src.cell, M)
        for ph = 1:nph
            q_ph = phase_source(c, src, rhoS[ph], kr, mu, ph)
            @inbounds acc[ph, src.cell] -= q_ph
        end
    end
end

# function insert_phase_sources!(acc::CuArray, model, kr, mu, rhoS, sources)
#     sources::CuArray
#     ix = map(cell, sources)
#     if !isa(rhoS, CuArray)
#         @warn "SurfaceDensities is not a CuArray, will convert whenever needed. Improve performance by converting once."  maxlog=1
#         rhoS = CuArray(rhoS)
#     end
#     rhoS::CuArray
#     @tullio acc[ph, ix[i]] = acc[ph, ix[i]] - phase_source(sources[i].cell, sources[i], rhoS[ph], kr, mu, ph)
# end

function convergence_criterion(model::SimulationModel{D, S}, storage, eq::ConservationLaw{:TotalMasses}, eq_s, r; dt = 1, update_report = missing) where {D, S<:MultiPhaseSystem}
    M = global_map(model.domain)
    v = x -> as_value(Jutul.active_view(x, M, for_variables = false))
    Φ = v(storage.state.FluidVolume)
    ρ = v(storage.state.PhaseMassDensities)

    nph = number_of_phases(model.system)
    cnv, mb = cnv_mb_errors(r, Φ, ρ, dt, Val(nph))

    names = phase_names(model.system)
    R = (CNV = (errors = cnv, names = names),
         MB = (errors = mb, names = names))
    return R
end

function cnv_mb_errors(r, Φ, ρ, dt, ::Val{N}) where N
    nc = length(Φ)
    mb = @MVector zeros(N)
    cnv = @MVector zeros(N)
    avg_density = @MVector zeros(N)

    pv_t = 0.0
    @inbounds for c in 1:nc
        pv_c = Φ[c]
        pv_t += pv_c
        @inbounds for ph = 1:N
            r_ph = r[ph, c]
            ρ_ph = ρ[ph, c]
            # MB
            mb[ph] += r_ph
            avg_density[ph] += abs(ρ_ph)
            # CNV
            cnv[ph] = max(cnv[ph], dt*abs(r_ph)/(ρ_ph*pv_c))
        end
    end
    @inbounds for ph = 1:N
        ρ_avg = avg_density[ph]/nc
        mb[ph] = (dt/pv_t)*abs(mb[ph])/ρ_avg
    end
    return (Tuple(cnv), Tuple(mb))
end

function cpr_weights_no_partials!(w, model::SimulationModel{R, S}, state, r, n, bz, scaling) where {R, S<:ImmiscibleSystem}
    ρ = state.PhaseMassDensities
    nc = size(w, 2)
    tb = minbatch(model.context, nc)
    M = global_map(model.domain)
    density = Jutul.active_view(ρ, M, for_variables = false)
    @batch minbatch = tb for i in 1:n
        for ph in axes(w, 1)
            @inbounds w[ph, i] = 1/value(density[ph, i])
        end
    end
end

function capillary_pressure(model, s)
    pck = :CapillaryPressure
    ref_index = get_reference_phase_index(model.system)
    if haskey(s, pck)
        pc = s[pck]
    else
        pc = nothing
    end
    return(pc, ref_index)
end

function get_reference_phase_index(::SinglePhaseSystem)
    return 1
end

"""
    get_reference_phase_index(system::JutulSystem)

Get the index of the reference phase in the system. The reference phase is the
pressure for which the pressure is given in the system. For single-phase systems
and models without capillary pressure this is unambigious.
"""
function get_reference_phase_index(system::JutulSystem)
    mphases = get_phases(system)
    return get_reference_phase_index(mphases)
end

function get_reference_phase_index(sys::MultiPhaseSystem)
    return sys.reference_phase_index
end

"""
    get_reference_phase_index(mphases)

Get the index of the reference phase for a set of phases. The reference phase is
selected as the liquid phase if present, otherwise the aqueous phase. If neither
is present the first phase is selected.
"""
function get_reference_phase_index(mphases)
    function find_phase(k)
        ix = 0
        for (i, ph) in enumerate(mphases)
            if ph isa k
                ix = i
                break
            end
        end
        return ix
    end
    l = find_phase(LiquidPhase)
    a = find_phase(AqueousPhase)
    if l > 0
        phase = l
    elseif a > 0
        phase = a
    else
        # Last phase is water or something custom, send out 1.
        phase = 1
    end
    return phase
end
