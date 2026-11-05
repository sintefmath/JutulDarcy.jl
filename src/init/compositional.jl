"COMPVD composition, with continuous total composition for EQUIL mode 1."
struct COMPVDComposition{V, L, D, F}
    vapor::V
    liquid::L
    contact::D
    overall::F
    mode::Int
end

function (composition::COMPVDComposition)(depth, phase::AbstractPhase)
    composition.mode == 1 && return composition.overall(depth)
    curve = phase isa VaporPhase ? composition.vapor : composition.liquid
    if isnothing(curve)
        curve = phase isa VaporPhase ? composition.liquid : composition.vapor
    end
    return curve(depth)
end

function (composition::COMPVDComposition)(depth)
    composition.mode == 1 && return composition.overall(depth)
    phase = depth < composition.contact ? VaporPhase() : LiquidPhase()
    return composition(depth, phase)
end

function compvd_composition_function(table, contact; normalize = :sum, mode = 1)
    n = size(table, 2) - 3
    flags = table[:, end-1]
    all(x -> x == 0 || x == 1, flags) || throw(ArgumentError("COMPVD phase flags must be 0 (vapor) or 1 (liquid)."))
    function curve(rows)
        isempty(rows) && return nothing
        depths = table[rows, 1]
        all(diff(depths) .> 0) || throw(ArgumentError("COMPVD depths must increase."))
        fractions = map(rows) do row
            z = SVector{n, Float64}(table[row, 2:n+1])
            all(isfinite, z) || throw(ArgumentError("COMPVD compositions must be finite."))
            if normalize == :clampsum
                z = clamp.(z, 0.0, 1.0)
            end
            all(x -> x >= 0, z) && sum(z) > 0 || throw(ArgumentError("COMPVD compositions must be nonnegative with a positive sum."))
            if normalize == :sum || normalize == :clampsum || normalize == true
                z /= sum(z)
            elseif normalize == false
                sum(z) ≈ 1 || throw(ArgumentError("COMPVD compositions must sum to one."))
            else
                throw(ArgumentError("Unknown composition normalization $normalize."))
            end
            return z
        end
        # The saturation-pressure column is observed data, not an imposed
        # pressure: hydrostatic pressures still come from EQUIL and the EOS.
        return get_1d_interpolator(depths, fractions)
    end
    return COMPVDComposition(curve(findall(isequal(0), flags)),
        curve(findall(isequal(1), flags)), contact, curve(axes(table, 1)), Int(mode))
end

function li_critical_temperature(eos, z)
    props = eos.mixture.properties
    volume = sum(i -> props[i].V_c*z[i], eachindex(z))
    return sum(i -> props[i].T_c*props[i].V_c*z[i], eachindex(z))/volume
end

function compositional_contact_state(eos, T, z, mode)
    MCF = MultiComponentFlash
    flash(p) = MCF.flashed_mixture_2ph(eos, (p=p, T=T, z=Vector{Float64}(z)); tolerance=1e-10, maxiter=300)
    upper = 1e8
    previous = upper
    lower = NaN
    for p in exp.(range(log(upper), log(1e4), length=120))
        f = flash(p)
        if f.state == MCF.two_phase_lv
            lower, upper = p, previous
            break
        end
        previous = p
    end
    isfinite(lower) || throw(ArgumentError("No saturation pressure found for EQUIL option $mode at T=$T."))
    for _ in 1:70
        p = (upper+lower)/2
        if flash(p).state == MCF.two_phase_lv
            lower = p
        else
            upper = p
        end
        upper-lower < 1e-9*upper && break
    end
    f = flash(lower*(1-1e-7))
    expected = mode == 3 ? f.V < 0.5 : f.V > 0.5
    expected || throw(ArgumentError("EQUIL option $mode specifies the wrong phase composition at the gas-oil contact."))
    return (upper+lower)/2, f.liquid.mole_fractions, f.vapor.mole_fractions
end
