module AFISetup
    using Jutul, JutulDarcy, GeoEnergyIO, OrderedCollections
    import GeoEnergyIO.IXParser: read_afi_file, AFIInputFile, find_records

    include("afi_reservoir.jl")
    include("afi_wells.jl")
    include("afi_system.jl")
    include("afi_init.jl")
    include("afi_pvt.jl")
    include("afi_saturation.jl")
end
