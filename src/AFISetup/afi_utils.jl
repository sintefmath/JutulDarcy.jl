
function merge_records(vals)
    I = firstindex(vals)
    kw = vals[I].keyword
    a = vals[I].value
    for i in eachindex(vals)
        if i == I
            continue
        end
        @assert vals[i].keyword == kw "Cannot merge dissimilar keywords: '$kw' and '$(vals[i].keyword)'"
        b = vals[i].value
        a = merge_records(a, b)
    end
    return (value = a, keyword = kw)
end

function merge_records(a, b)
    keysa = keys(a)
    keysb = keys(b)
    out = Dict{String, Any}()
    for k in intersect(keysa, keysb)
        out[k] = merge(a[k], b[k])
    end
    for ka in setdiff(keysa, keysb)
        out[ka] = a[ka]
    end
    for kb in setdiff(keysb, keysa)
        out[kb] = b[kb]
    end
    return out
end
