export preprocess_relperm_table

function preprocess_relperm_table(swof, ϵ = 1e-16; swcon = 0.0)
    sw = vec(swof[:, 1])
    krw = vec(swof[:, 2])
    so = 1 .- sw
    so = vec(so[end:-1:1])
    if swcon != 0
        for i in eachindex(so)
            so[i] -= swcon
        end
    end
    kro = vec(swof[end:-1:1, 3])
    # Make sure that we don't extrapolate
    sw, krw = add_missing_endpoints(sw, krw)
    # Change so table to be with respect to so,
    # and to be increasing with respect to input
    so, kro = add_missing_endpoints(so, kro)
    # Subtract a tiny bit from the saturations at endpoints.
    # This is to ensure that the derivative ends up as zero
    # when evaluated at s corresponding to kr_max
    ensure_endpoints!(so, kro, ϵ)
    ensure_endpoints!(sw, krw, ϵ)
    s = [sw, so]
    krt = [krw, kro]
    return s, krt
end

function add_missing_endpoints(s, kr)
    if s[1] > 0.0
        s = vcat(0.0, s)
        kr = vcat(0.0, kr)
    end
    if s[end] < 1.0
        s = vcat(s, 1.0)
        kr = vcat(kr, kr[end])
    end
    return (s, kr)
end

function ensure_endpoints!(x, f, ϵ)
    n = length(x)
    for i in (n-1):-1:2
        if f[i] != f[i-1]
            x[i] -= ϵ
            break
        end
    end
    for i in 1:(n-1)
        if f[i] != f[i+1]
            x[i] -= ϵ
            break
        end
    end
end