module Geothermal

    using Jutul, JutulDarcy
    export setup_btes_well, setup_vertical_btes_well
    export BTESWellSupplyToReturnMassCT, ClosedLoopSupplyToReturnEnergyCT, BTESWellGroutEnergyCT
    export update_cross_term_in_entity!

    function JutulDarcy.setup_reservoir_model(reservoir::DataDomain, ::Val{:geothermal}; kwarg...)
        return setup_reservoir_model_geothermal(reservoir; kwarg...)
    end

    include("wells/btes.jl")
    include("wells/cross_terms.jl")
    include("properties.jl")

end
