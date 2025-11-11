module AFISetup
    using Jutul, JutulDarcy, GeoEnergyIO, OrderedCollections
    using Dates
    import GeoEnergyIO.IXParser: read_afi_file, AFIInputFile, find_records

    include("afi_reservoir.jl")
    include("afi_wells.jl")
    include("afi_system.jl")
    include("afi_init.jl")
    include("afi_pvt.jl")
    include("afi_saturation.jl")
    include("afi_regions.jl")
    include("afi_model.jl")
    include("afi_schedule.jl")
    include("afi_utils.jl")

    function setup_case_from_afi(x::AbstractDict; kwarg...)
        return setup_case_from_afi(AFIInputFile(x); kwarg...)
    end

    function setup_case_from_afi(afi::AFIInputFile; kwarg...)
        model, prm = setup_reservoir_model(afi; extra_out = true, kwarg...)
        state0 = setup_reservoir_state(afi, model)
        dt, forces = setup_afi_schedule(afi, model)
        return Jutul.JutulCase(model, dt, forces, state0 = state0, parameters = prm)
    end
end
