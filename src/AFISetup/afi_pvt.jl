function setup_pvt_variables(d::AFIInputFile, sys::Union{StandardBlackOilSystem, ImmiscibleSystem}, reservoir)
    bo_model = find_records(d, "BlackOilFluidModel", "IX", steps = false, model = true, once = false)
    if length(bo_model) == 0
        error("No BlackOilFluidModel record found in AFI file. This function is currently limited to blackoil type models.")
    end
    if length(bo_model) > 1
        @warn "Multiple BlackOilFluidModel records found in AFI file. Using the first one."
    end
    bo_model = bo_model[1].value
    phases = JutulDarcy.get_phases(sys)
    pvt = []
    if AqueousPhase() in phases
        wc = bo_model["WaterCompressibilities"]
        water_tab_raw = ConstMuBTable(wc["RefPressure"], 1.0/wc["FormationVolumeFactor"], wc["Compressibility"], wc["Viscosity"], wc["ViscosityCompressibility"])
        water_tab = JutulDarcy.PVTW((water_tab_raw, ))
        push!(pvt, water_tab)
    end
    if LiquidPhase() in phases
        if JutulDarcy.has_disgas(sys)
            dtab = bo_model["OilTable"]["table"]
            pvto_like = JutulDarcy.AFISetup.to_processed_pvt_table(dtab, "SolutionGOR", ["Pressure", "FormationVolumeFactor", "Viscosity"])
            pvto_ext = GeoEnergyIO.InputParser.restructure_pvt_table(pvto_like)
            oil_tab_raw = JutulDarcy.PVTOTable(pvto_ext)
            oil_tab = JutulDarcy.PVTO(oil_tab_raw)
        else
            dtab = bo_model["DeadOilTable"]["table"]
            p = dtab["Pressure"]
            Bo = dtab["FormationVolumeFactor"]
            mu = dtab["Viscosity"]
            oil_tab_raw = JutulDarcy.MuBTable(p, 1 ./ Bo, mu)
            oil_tab = JutulDarcy.PVCDO((oil_tab_raw, ))
        end
        push!(pvt, oil_tab)
    end
    if VaporPhase() in phases
        if JutulDarcy.has_vapoil(sys)
            @info "Not finished" bo_model["GasTable"]
            error("Not yet implemented")
        else
            dtab = bo_model["UndersaturatedGasTable"]["table"]
            p = dtab["Pressure"]
            Bg = dtab["FormationVolumeFactor"]
            mu = dtab["Viscosity"]
            gas_tab_raw = JutulDarcy.MuBTable(p, 1 ./ Bg, mu)
            gas_tab = JutulDarcy.PVDG((gas_tab_raw, ))
        end
        push!(pvt, gas_tab)
    end
    pvt_vars = Dict()
    pvt_vars[:PhaseMassDensities] = DeckPhaseMassDensities(pvt)
    pvt_vars[:ShrinkageFactors] = DeckShrinkageFactors(pvt)
    pvt_vars[:PhaseViscosities] = DeckPhaseViscosities(pvt)
    return pvt_vars
end

function rs_table_from_oil_table(x)
    tab = x["table"]
    ix = tab["SubTableIndex"]
    gor = tab["SolutionGOR"]
    p = tab["Pressure"]
    p_t, rs_t = extract_saturated_table(ix, p, gor)
    return JutulDarcy.saturated_table(p_t, rs_t)
end

function extract_saturated_table(table_indices, p, ratio)
    if eltype(table_indices) != Int
        table_indices = Int.(table_indices)
    end
    ratio_t = Float64[]
    p_t = Float64[]
    prev_table = -1
    for (row, table_index) in enumerate(table_indices)
        if table_index != prev_table
            push!(p_t, p[row])
            push!(ratio_t, ratio[row])
            prev_table = table_index
        end
    end
    return (p_t, ratio_t)
end

function rv_table_from_gas_table(x)
    tab = x["table"]
    ix = tab["SubTableIndex"]
    rv = tab["VaporOGR"]
    p = tab["Pressure"]
    p_t, rv_t =  extract_saturated_table(ix, p, rv)
    return JutulDarcy.saturated_table(p_t, rv_t)
end

function to_processed_pvt_table(dtab, main_key, key_order)
    tabix = Int.(dtab["SubTableIndex"])
    mintab = minimum(tabix)
    if mintab == 0
        @. tabix += 1
    else
        mintab == 1 || error("Unexpected SubTableIndex starting at $mintab")
    end
    issorted(tabix) || error("SubTableIndex not sorted")
    maxtab = maximum(tabix)
    out = Vector{Float64}[]
    for subtab_no in 1:maxtab
        next = Float64[]
        is_first = true
        for (i, ix) in enumerate(tabix)
            if ix != subtab_no
                continue
            end
            if is_first
                push!(next, dtab[main_key][i])
                is_first = false
            end
            for k in key_order
                push!(next, dtab[k][i])
            end
        end
        if is_first
            error("No entries found for table index $i")
        end
        push!(out, next)
    end
    return out
end

