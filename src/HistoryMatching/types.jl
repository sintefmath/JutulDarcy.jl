const I_type = Jutul.LinearInterpolant{Vector{Float64}, Vector{Float64}, Missing}

struct MatchPeriod
    start::Float64
    stop::Float64
    weights::Vector{Float64}
    step_idx::Vector{Int}
end

function MatchPeriod(start::Real, stop::Real)
    return MatchPeriod(start, stop, Float64[], Int[])
end

struct WellMatch
    name::Symbol
    type::String
    weight::Union{Float64, Vector{Float64}}
    scale::Float64
    is_injector::Bool
    data::I_type
end

struct ReservoirMatch
    name::Symbol
    weight::Float64
    scale::Float64
    row::Int
    data::I_type
end

mutable struct HistoryMatchLogger
    data::Union{Missing, Dict{Any, Any}}
end

function HistoryMatchLogger()
    return HistoryMatchLogger(missing)
end

struct HistoryMatch
    case::JutulCase
    states
    summary
    periods::Union{Missing, Vector{MatchPeriod}}
    wellpos::Dict{Symbol, Int}
    injector_rate::Vector{WellMatch}
    injector_orat::Vector{WellMatch}
    injector_wrat::Vector{WellMatch}
    injector_grat::Vector{WellMatch}
    injector_bhp::Vector{WellMatch}
    producer_rate::Vector{WellMatch}
    producer_bhp::Vector{WellMatch}
    producer_grat::Vector{WellMatch}
    producer_orat::Vector{WellMatch}
    producer_lrat::Vector{WellMatch}
    producer_wrat::Vector{WellMatch}
    producer_cumulative_oil::Vector{WellMatch}
    producer_cumulative_gas::Vector{WellMatch}
    producer_cumulative_water::Vector{WellMatch}
    producer_cumulative_liquid::Vector{WellMatch}
    reservoir::Vector{ReservoirMatch}
    total_scale::Float64
    logger::HistoryMatchLogger
end

function HistoryMatch(case::JutulCase; kwarg...)
    return HistoryMatch(case, missing, missing; kwarg...)
end

function HistoryMatch(case::JutulCase, res::JutulDarcy.ReservoirSimResult, smry = missing; kwarg...)
    if ismissing(smry)
        smry = JutulDarcy.summary_result(case, res, :si, field = false, wells = true)
    end
    states = res.states
    return HistoryMatch(case, states, smry; kwarg...)
end

function HistoryMatch(case::JutulCase, states, summary;
        periods = missing,
        period_weights = missing,
        normalized_periods = missing,
        scale = 1.0#/sum(case.dt)
    )
    periods = setup_periods(case, periods, period_weights, normalized_periods)
    case.model::MultiModel
    wellpos = Dict{Symbol, Int}()
    for w in get_wells(case)
        wellpos[w] = JutulDarcy.get_well_position(case.model[:Facility].domain, w)
    end
    return HistoryMatch(
        case,
        states,
        summary,
        periods,
        wellpos,
        WellMatch[],  # injector_rate
        WellMatch[],  # injector_orat
        WellMatch[],  # injector_wrat
        WellMatch[],  # injector_grat
        WellMatch[],  # injector_bhp
        WellMatch[],  # producer_rate
        WellMatch[],  # producer_bhp
        WellMatch[],  # producer_grat
        WellMatch[],  # producer_orat
        WellMatch[],  # producer_lrat
        WellMatch[],  # producer_wrat
        WellMatch[],  # producer_cumulative_oil
        WellMatch[],  # producer_cumulative_gas
        WellMatch[],  # producer_cumulative_water
        WellMatch[],  # producer_cumulative_liquid
        ReservoirMatch[],  # reservoir
        scale,
        HistoryMatchLogger()
    )
end

struct GlobalHistoryMatchObjective <: Jutul.AbstractGlobalObjective
    match::HistoryMatch
    evaluation_count
end

function Base.show(io::IO, obj::GlobalHistoryMatchObjective)
    print(io, "GlobalHistoryMatchObjective:\n")
    Base.show(io, obj.match)
end

struct SumHistoryMatchObjective <: Jutul.AbstractSumObjective
    match::HistoryMatch
    evaluation_count
end

function Base.show(io::IO, obj::SumHistoryMatchObjective)
    println(io, "SumHistoryMatchObjective")
    Base.show(io, obj.match)
end

function Base.show(io::IO, hm::HistoryMatch)
    tstr = Jutul.get_tstr(sum(hm.case.dt))
    print(io, "HistoryMatch objective covering case with $tstr total simulation time\n")
    if !get(io, :compact, false)
        println("")
        keys = propertynames(hm)
        # Match, well, injector/producer, weight, scale
        wdata = Vector{WellMatch}()
        for k in keys
            v = getfield(hm, k)
            if v isa Vector{WellMatch}
                append!(wdata, v)
            end
        end
        sort!(wdata, by = x -> "$(x.is_injector ? "I" : "P")_$(x.name)_$(x.type)")

        header = ["Well", "Type", "Quantity", "Weight", "Scale"]
        tab = Matrix{Any}(undef, length(wdata), length(header))
        for (i, v) in enumerate(wdata)
            if v.is_injector
                nm = "Injector"
            else
                nm = "Producer"
            end
            tab[i, 1] = v.name
            tab[i, 2] = nm
            tab[i, 3] = v.type
            tab[i, 4] = v.weight
            tab[i, 5] = v.scale
        end
        pretty_table(io, tab; header = header, crop = :horizontal)
    end
end


const HistoryMatchObjective = Union{GlobalHistoryMatchObjective, SumHistoryMatchObjective}

