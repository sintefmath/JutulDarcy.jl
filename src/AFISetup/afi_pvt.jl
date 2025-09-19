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
        water_tab = ConstMuBTable(wc["RefPressure"], 1.0/wc["FormationVolumeFactor"], wc["Compressibility"], wc["Viscosity"], wc["ViscosityCompressibility"])
        push!(pvt, water_tab)
    end
    if LiquidPhase() in phases
        if JutulDarcy.has_disgas(sys)
            dtab = bo_model["OilTable"]["table"]
            p_bub = sys.rs_max[1].X
            pos = Int.(dtab["SubTableIndex"]) .+ 1
            rs = dtab["SolutionGOR"]
            pressure = dtab["Pressure"]
            shrinkage = 1 ./ dtab["FormationVolumeFactor"]
            viscosity = dtab["Viscosity"]
            oil_tab = JutulDarcy.PVTOTable(pos, rs, pressure, p_bub, shrinkage, viscosity)
        else
            dtab = bo_model["DeadOilTable"]["table"]
            p = dtab["Pressure"]
            Bo = dtab["FormationVolumeFactor"]
            mu = dtab["Viscosity"]
            oil_tab = JutulDarcy.MuBTable(p, 1 ./ Bo, mu)
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
            gas_tab = JutulDarcy.MuBTable(p, 1 ./ Bg, mu)
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
    p_t = Float64[]
    rs_t = Float64[]
    prev = 0
    for idx in ix
        idx = Int(idx)
        if idx != prev
            push!(p_t, p[idx])
            push!(rs_t, gor[idx])
            prev = idx
        end
    end
    return JutulDarcy.saturated_table(p_t, rs_t)
end

function rv_table_from_gas_table(x)
    tab = x["table"]
    ix = tab["SubTableIndex"]
    rv = tab["VaporOGR"]
    p = tab["Pressure"]
    p_t = Float64[]
    rv_t = Float64[]
    prev = 0
    for idx in ix
        idx = Int(idx)
        if idx != prev
            push!(p_t, p[idx])
            push!(rv_t, rv[idx])
            prev = idx
        end
    end
    return JutulDarcy.saturated_table(p_t, rs_t)
end
