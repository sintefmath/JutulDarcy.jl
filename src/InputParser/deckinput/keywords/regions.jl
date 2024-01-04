function finish_current_section!(data, units, cfg, outer_data, ::Val{:REGIONS})
    
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:RPTREGS})
    parser_message(cfg, outer_data, "RPTREGS", PARSER_MISSING_SUPPORT)
    rec = read_record(f)
    # TODO: Process the record.
end
