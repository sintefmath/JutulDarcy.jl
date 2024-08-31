module CO2Properties
    using Jutul, JutulDarcy, StaticArrays, MultiComponentFlash, DelimitedFiles, Artifacts, LazyArtifacts
    include("kvalues.jl")
    include("props.jl")
    include("reading.jl")
    include("setup.jl")

    function JutulDarcy.setup_reservoir_model(reservoir::DataDomain, ::Val{:co2brine}; kwarg...)
        return setup_reservoir_model_co2_brine(reservoir; kwarg...)
    end
end
