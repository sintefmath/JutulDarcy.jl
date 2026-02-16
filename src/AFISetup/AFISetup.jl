module AFISetup
    using Jutul, JutulDarcy, GeoEnergyIO, OrderedCollections
    using Dates
    import GeoEnergyIO.IXParser: read_afi_file, AFIInputFile, find_records

    export setup_case_from_afi

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
    include("interface.jl")
end
