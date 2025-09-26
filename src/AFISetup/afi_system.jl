
function setup_system(d::AFIInputFile)
    phases = setup_phases(d)
    bo_model = find_records(d, "BlackOilFluidModel", "IX", steps = false, model = true, once = true)
    has_water = AqueousPhase() in phases
    has_oil = LiquidPhase() in phases
    has_gas = VaporPhase() in phases

    if isnothing(bo_model)
        error("No BlackOilFluidModel record found in AFI file. This function is currently limited to blackoil type models.")
    else
        bo_model = bo_model.value
        has_disgas = haskey(bo_model, "OilTable")
        has_vapoil = haskey(bo_model, "GasTable")
        has_vapoil && error("Vaporizing oil models are not yet supported in setup from AFI.")
        rs_max = rv_max = nothing
        if has_disgas
            rs_max = rs_table_from_oil_table(bo_model["OilTable"])
        end
        if has_vapoil
            rv_max = rv_table_from_gas_table(bo_model["GasTable"])
        end
        rhoS = Float64[]
        if has_water
            push!(rhoS, bo_model["WaterCompressibilities"]["SurfaceDensity"])
        end
        if has_oil
            push!(rhoS, bo_model["OilSurfaceDensity"])
        end
        if has_gas
            push!(rhoS, bo_model["GasSurfaceDensity"])
        end
        if has_disgas || has_vapoil
            sys = StandardBlackOilSystem(
                phases = phases,
                rs_max = rs_max,
                rv_max = rv_max,
                reference_densities = rhoS
            )
        else
            sys = ImmiscibleSystem(phases, reference_densities = rhoS)
        end
    end
    return sys
end

function setup_phases(d::AFIInputFile)
    phases = find_records(d, "PhasesPresent", "IX", steps = false, model = true, once = true)
    if isnothing(phases)
        sim = find_records(d, "Simulation", "IX", steps = false, model = true, once = true)
        if !isnothing(sim)
            phases = find_records(sim.value.value, "PhasesPresent", once = true)
            if !isnothing(phases)
                phases = map(String, phases.value)
            end
        end
    else
        phases = phases.value
    end
    if isnothing(phases)
        println("Did not find phases - assuming water, oil and gas")
        phases = ["WATER", "OIL", "GAS"]
    end

    phases = map(x -> strip(uppercase(x)), phases)
    length(setdiff(phases, ["GAS", "OIL", "WATER"])) == 0 || error("Only GAS, OIL, and WATER phases are supported in AFI files.")
    # Make sure that the phases are in the canonical order for JutulDarcy (WOG)
    # In theory this is not needed, but it will simplify well controls.
    typed_phases = []
    if "WATER" in phases
        push!(typed_phases, AqueousPhase())
    end
    if "OIL" in phases
        push!(typed_phases, LiquidPhase())
    end
    if "GAS" in phases
        push!(typed_phases, VaporPhase())
    end
    return Tuple(typed_phases)
end
