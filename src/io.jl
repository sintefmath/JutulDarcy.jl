export table_to_relperm

function table_to_relperm(swof; swcon = 0.0, first_label = :w, second_label = :ow)
    sw = vec(swof[:, 1])
    krw = vec(swof[:, 2])
    krw = PhaseRelPerm(sw, krw, label = first_label)
    kro = vec(swof[end:-1:1, 3])
    so = 1 .- sw
    so = vec(so[end:-1:1])
    @. so = so - swcon
    krow = PhaseRelPerm(so, kro, label = second_label)
    return (krw, krow)
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