
function parse_keyword!(data, outer_data, units, f, ::Val{:NSTACK})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:RPTSOL})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:EXCEL})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:RUNSUM})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:SEPARATE})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:NEWTON})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:TCPU})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:ELAPSED})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:MAXDPR})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:FOPR})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:FWPR})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:FWIR})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:FOPT})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:FWPT})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:FWIT})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:RPTONLY})
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WMCTL})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:RPTRST})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WLPR})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WGOR})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WWIR})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WWPR})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WOPR})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WBHP})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:BGSAT})
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:BOSAT})
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:BOKR})
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:BWSAT})
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:BPR})
    skip_record(f)
end
