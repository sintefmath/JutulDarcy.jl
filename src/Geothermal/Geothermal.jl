module Geothermal

    using Jutul, JutulDarcy
    export setup_btes_well, setup_vertical_btes_well
    export BTESWellSupplyToReturnMassCT, BTESWellSupplyToReturnEnergyCT, BTESWellGroutEnergyCT
    export update_cross_term_in_entity!

    include("wells/btes.jl")
    include("wells/cross_terms.jl")

end