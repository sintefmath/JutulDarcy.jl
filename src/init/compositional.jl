"Composition curves separated at the gas-oil contact specified by EQUIL."
struct COMPVDComposition{V, L, D}
    vapor::V
    liquid::L
    contact::D
end

function (composition::COMPVDComposition)(depth, phase::AbstractPhase)
    curve = phase isa VaporPhase ? composition.vapor : composition.liquid
    if isnothing(curve)
        curve = phase isa VaporPhase ? composition.liquid : composition.vapor
    end
    return curve(depth)
end

function (composition::COMPVDComposition)(depth)
    phase = depth < composition.contact ? VaporPhase() : LiquidPhase()
    return composition(depth, phase)
end

function compvd_composition_function(table, contact; normalize = :sum)
    n = size(table, 2) - 3
    flags = table[:, end-1]
    all(x -> x == 0 || x == 1, flags) || throw(ArgumentError("COMPVD phase flags must be 0 (vapor) or 1 (liquid)."))
    function phase_curve(flag)
        rows = findall(isequal(flag), flags)
        isempty(rows) && return nothing
        depths = table[rows, 1]
        all(diff(depths) .> 0) || throw(ArgumentError("COMPVD depths must increase within each phase."))
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
    return COMPVDComposition(phase_curve(0), phase_curve(1), contact)
end
