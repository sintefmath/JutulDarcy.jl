module AFISetup
    using Jutul, JutulDarcy, GeoEnergyIO, OrderedCollections
    import GeoEnergyIO.IXParser: read_afi_file, AFIInputFile, find_records

    include("reservoir.jl")
    include("wells.jl")
end
