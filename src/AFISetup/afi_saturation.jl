function setup_saturation_variables(d::AFIInputFile, reservoir, sys)
    phases = JutulDarcy.get_phases(sys)
    nph = length(phases)
    has_water = AqueousPhase() in phases
    has_oil = LiquidPhase() in phases
    has_gas = VaporPhase() in phases

    satfuns = find_records(d, "SaturationFunction", "IX", steps = false, model = true)
end
