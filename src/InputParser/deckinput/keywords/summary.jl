function finish_current_section!(data, units, cfg, outer_data, ::Val{:SUMMARY})
    
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:NSTACK})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:RPTSOL})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:EXCEL})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:RUNSUM})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:SEPARATE})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:NEWTON})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:TCPU})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:ELAPSED})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MAXDPR})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FOPR})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FWPR})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FWIR})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FOPT})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FWPT})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FWIT})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:RPTONLY})
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WMCTL})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:RPTRST})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WLPR})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WGOR})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WWIR})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WWPR})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WOPR})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WBHP})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:BGSAT})
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:BOSAT})
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:BOKR})
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:BWSAT})
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:BPR})
    skip_record(f)
end
