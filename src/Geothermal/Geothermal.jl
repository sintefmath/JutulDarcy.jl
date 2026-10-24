module Geothermal

    using Jutul, JutulDarcy

    function JutulDarcy.setup_reservoir_model(reservoir::DataDomain, ::Val{:geothermal}; kwarg...)
        return setup_reservoir_model_geothermal(reservoir; kwarg...)
    end

    include("properties.jl")

end
