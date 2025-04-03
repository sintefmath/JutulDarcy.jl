"""
    compute_contraction_factor(r, N)

Compute contraction factor from a number of Newton iterate distances from
convergence (defined with a user-prescribed distance metric). For more than two
iterates, the contraction factor is estimated using least-squares assuming the
iterates follow a geometric series: rₙ = Θⁿr₀. The function also computes the
target contraction factor under this assumption assuming N iterations left.
"""
function compute_contraction_factor(d, N)

    # Use distance from convergence + 1 to avoid division by small numbers
    r = d .+ 1.0
    # Number of iterates used to estimate contraction factor
    n = size(r,1)-1
    # Approximate log10 of contraction factor using least squares
    num = sum(log.(r[2:end]./r[1]).*(1:n))
    den = sum((1:n).^2)
    θ = exp(num./den)
    # Compute target contraction factor convergence in N iterations
    θ_target = r[1].^(-1/N)
    
    return θ, θ_target

end

function oscillation(contraction_factors, tol)

    θ = contraction_factors
    (length(θ) < 3) ? (return false) : nothing

    θ_1 = θ[end-2]
    θ_2 = θ[end-1]
    θ_3 = θ[end]

    ok_1 = θ_1 .< tol
    ok_2 = θ_2 .< tol
    ok_3 = θ_3 .< tol

    return .&(xor.(ok_1, ok_2), xor.(ok_2, ok_3))

end