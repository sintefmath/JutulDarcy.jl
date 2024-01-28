function current_unit_system(deck)
    rs = deck["RUNSPEC"]
    choices = ["METRIC", "SI", "LAB", "FIELD"]
    found = false
    sys = missing
    for choice in choices
        if haskey(rs, choice)
            @assert ismissing(sys) "Cannot have multiple unit keywords (encountered $choice and $sys)"
            sys = choice
        end
    end
    if ismissing(sys)
        @warn "Units not found. Assuming field units."
        out = :field
    else
        out = Symbol(lowercase(sys))
    end
    return out
end

Base.@kwdef struct DeckUnitSystem{S, T}
    length::T = 1.0
    area::T = 1.0
    time::T = 1.0
    density::T = 1.0
    pressure::T = 1.0
    mol::T = 1.0
    mass::T = 1.0
    u_rs::T = 1.0
    u_rv::T = 1.0
    concentration::T = 1.0
    compressibility::T = 1.0
    viscosity::T = 1.0
    surface_tension::T = 1.0
    jsurface_tension::T = 1.0
    permeability::T = 1.0
    liquid_volume_surface::T = 1.0
    liquid_volume_reservoir::T = 1.0
    liquid_formation_volume_factor::T = 1.0
    gas_volume_surface::T = 1.0
    gas_volume_reservoir::T = 1.0
    gas_formation_volume_factor::T = 1.0
    volume::T = 1.0
    transmissibility::T = 1.0
    rock_conductivity::T = 1.0
    volume_heat_capacity::T = 1.0
    mass_heat_capacity::T = 1.0
    molar_mass::T = 1.0
    relative_temperature::Symbol = :Celsius
    absolute_temperature::Symbol = :Kelvin
end

function DeckUnitSystem(sys::Symbol, T = Float64)
    #meter, day, kilogram, bar = si_units(:meter, :day, :kilogram, :bar)
    u = Jutul.all_units()
    m = u[:meter]
    K = u[:kelvin]
    day = u[:day]
    centi = u[:centi]
    kilogram = u[:kilogram]
    ft = u[:feet]
    psi = u[:psi]
    pound = u[:pound]
    kilo = u[:kilo]
    stb = u[:stb]
    rankine = u[:rankine]
    btu = u[:btu]

    # Commons
    cP = u[:centi]*u[:poise]
    mD = u[:milli]*u[:darcy]
    if sys == :metric
        len = m
        kJ = u[:kilo]*u[:joule]
        volume = m^3
        time = day
        pressure = u[:bar]
        mol = u[:kilo]
        molar_mass = 1/(mol)
        mass = kilogram
        viscosity = cP
        surface_tension = u[:newton]/m
        jsurface_tension = u[:dyne]/(centi*m)
        permeability = mD
        liquid_volume_surface = volume
        liquid_volume_reservoir = volume
        gas_volume_surface = volume
        gas_volume_reservoir = volume
        rock_conductivity = kJ/(m*day*K)
        volume_heat_capacity = kJ/(volume*K)
        mass_heat_capacity = kJ/(mass*K)
        relative_temperature = :Celsius
        absolute_temperature = :Kelvin
    elseif sys == :field
        len = ft
        time = day
        pressure = psi
        mol = pound*kilo
        molar_mass = 1/(mol)
        mass = pound
        viscosity = cP
        surface_tension = u[:lbf]/u[:inch]
        jsurface_tension = u[:dyne]/(centi*m)
        permeability = mD
        liquid_volume_surface = stb
        liquid_volume_reservoir = stb
        gas_volume_surface = kilo*ft^3
        gas_volume_reservoir = stb
        rock_conductivity = btu / (ft*day*rankine)
        volume_heat_capacity = btu / (ft^3*rankine)
        mass_heat_capacity = btu / (pound*rankine)
        relative_temperature = :Fahrenheit
        absolute_temperature = :Rankine
    elseif sys == :lab
        error("Not implemented")
    else
        @assert sys == :si
        len = 1.0
        time = 1.0
        pressure = 1.0
        mol = 1.0
        mass = 1.0
        viscosity = 1.0
        surface_tension = 1.0
        jsurface_tension = 1.0
        permeability = 1.0
        liquid_volume_surface = 1.0
        liquid_volume_reservoir = 1.0
        gas_volume_surface = 1.0
        gas_volume_reservoir = 1.0
        rock_conductivity = 1.0
        volume_heat_capacity = 1.0
        mass_heat_capacity = 1.0
        molar_mass = 1.0
        relative_temperature = :Celsius
        absolute_temperature = :Kelvin
    end
    area = len^2
    volume = len^3
    density = mass/volume
    concentration = mass/volume
    compressibility = 1.0/pressure
    transmissibility = viscosity * liquid_volume_reservoir / (time * pressure)
    return DeckUnitSystem{sys, T}(
        length = len,
        area = area,
        time = time,
        u_rs = gas_volume_surface/liquid_volume_surface,
        u_rv = liquid_volume_surface/gas_volume_surface,
        density = density,
        pressure = pressure,
        mol = mol,
        mass = mass,
        concentration = concentration,
        compressibility = compressibility,
        viscosity = viscosity,
        surface_tension = surface_tension,
        jsurface_tension = jsurface_tension,
        permeability = permeability,
        liquid_volume_surface = liquid_volume_surface,
        liquid_volume_reservoir = liquid_volume_reservoir,
        liquid_formation_volume_factor = liquid_volume_reservoir/liquid_volume_surface,
        gas_volume_surface = gas_volume_surface,
        gas_volume_reservoir = gas_volume_reservoir,
        gas_formation_volume_factor = gas_volume_reservoir/gas_volume_surface,
        volume = volume,
        molar_mass = molar_mass,
        transmissibility = transmissibility,
        rock_conductivity = rock_conductivity,
        volume_heat_capacity = volume_heat_capacity,
        mass_heat_capacity = mass_heat_capacity,
        relative_temperature = relative_temperature,
        absolute_temperature = absolute_temperature
    )
end

function deck_unit_system_label(::DeckUnitSystem{S, T}) where {S, T}
    return S
end

function swap_unit_system_axes!(x::AbstractMatrix, systems, eachunit; dim = 2)
    @assert eltype(eachunit)<:Symbol
    @assert size(x, dim) == length(eachunit)
    if dim == 1
        x_t = x'
    else
        x_t = x
    end
    for i in axes(x_t, 2)
        x_i = view(x_t, :, i)
        swap_unit_system!(x_i, systems, eachunit[i])
    end
    return x
end

function swap_unit_system_axes!(x::AbstractVector, systems, eachunit)
    @assert eltype(eachunit)<:Symbol
    @assert length(x) == length(eachunit) "Recieved vector of length $(length(x)) but units were $(length(eachunit)) long."
    for i in eachindex(x)
        x[i] = swap_unit_system(x[i], systems, eachunit[i])
    end
    return x
end

function swap_unit_system!(x::AbstractArray, systems, k)
    return swap_unit_system!(x, systems, Val(k))
end

function swap_unit_system!(x::AbstractArray, systems, k::Val)
    for i in eachindex(x)
        x[i] = swap_unit_system(x[i], systems, k)
    end
    return x
end

function swap_unit_system(val, systems, k::Symbol)
    return swap_unit_system(val, systems, Val(k))
end

function swap_unit_system(v, systems::Union{Nothing, Missing}, k::Val)
    # No systems - trivial conversion
    return v
end

function swap_unit_system(v, systems::Union{Nothing, Missing}, k::Symbol)
    # No systems - trivial conversion
    return v
end

function swap_unit_system(val, systems::NamedTuple, ::Union{Val{:identity}, Val{:id}})
    # Identity specifically means no unit.
    return val
end

function identity_unit_vector(x)
    return identity_unit_vector(length(x))
end

function identity_unit_vector(n::Int)
    utypes = Vector{Symbol}(undef, n)
    fill!(utypes, :id)
    return utypes
end

function swap_unit_system(val, systems::NamedTuple, ::Val{k}; reverse = false) where k
    (; to, from) = systems
    if reverse
        to, from = from, to
    end
    to_unit = deck_unit(to, k)
    from_unit = deck_unit(from, k)

    val_si = convert_to_si(val, from_unit)
    val_final = convert_from_si(val_si, to_unit)
    return val_final
end

function deck_unit(sys::DeckUnitSystem, s::Symbol)
    return deck_unit(sys, Val(s))
end

function deck_unit(sys::DeckUnitSystem, ::Val{k}) where k
    return getproperty(sys, k)
end

# Magic type overloads

function deck_unit(sys::DeckUnitSystem, ::Val{:Kh})
    return deck_unit(sys, :permeability)*deck_unit(sys, :length)
end

function deck_unit(sys::DeckUnitSystem, ::Val{:gigapascal})
    return si_unit(:Pa)*si_unit(:giga)
end

function deck_unit(sys::DeckUnitSystem, ::Val{:time_over_volume})
    return deck_unit(sys, :time)/deck_unit(sys, :volume)
end

function deck_unit(sys::DeckUnitSystem, ::Val{:liquid_rate_surface})
    return deck_unit(sys, :liquid_volume_surface)/deck_unit(sys, :time)
end

function deck_unit(sys::DeckUnitSystem, ::Val{:gas_rate_surface})
    return deck_unit(sys, :gas_volume_surface)/deck_unit(sys, :time)
end

function deck_unit(sys::DeckUnitSystem, ::Val{:liquid_rate_reservoir})
    return deck_unit(sys, :liquid_volume_reservoir)/deck_unit(sys, :time)
end

function deck_unit(sys::DeckUnitSystem, ::Val{:gas_rate_reservoir})
    return deck_unit(sys, :gas_volume_reservoir)/deck_unit(sys, :time)
end

function deck_unit(sys::DeckUnitSystem, ::Val{:critical_volume})
    return deck_unit(sys, :volume)/deck_unit(sys, :mol)
end
