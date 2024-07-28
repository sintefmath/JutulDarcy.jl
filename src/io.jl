function table_to_relperm(swof; swcon = 0.0, first_label = :w, second_label = :ow)
    sw = vec(swof[:, 1])
    krw = vec(swof[:, 2])
    krw = PhaseRelativePermeability(sw, krw, label = first_label)
    kro = vec(swof[end:-1:1, 3])
    so = 1 .- sw
    so = vec(so[end:-1:1])
    @. so = so - swcon
    krow = PhaseRelativePermeability(so, kro, label = second_label)
    return (krw, krow)
end

function saturation_table_handle_defaults(s, f)
    if any(isnan, f)
        # NaN values are removed due to INPUT file shenanigans
        ix = findall(!isnan, f)
        s = s[ix]
        f = f[ix]
    end
    return (s, f)
end

function add_missing_endpoints(s, kr)
    copied = false
    if s[1] > 0.0
        copied = true
        s = vcat(0.0, s)
        kr = vcat(0.0, kr)
    end
    if s[end] < 1.0
        copied = true
        s = vcat(s, 1.0)
        kr = vcat(kr, kr[end])
    end
    if !copied
        s = copy(s)
        kr = copy(kr)
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