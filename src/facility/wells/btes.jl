function setup_btes_well(D::DataDomain, reservoir_cells; name = :BTES, kwarg...)

    supply_well = setup_well(D::DataDomain, reservoir_cells;
        name = Symbol(name, "_supply"), WI = 0.0, simple_well = false, type = :btes, kwarg...)

    return_well = setup_well(D::DataDomain, reservoir_cells;
        name = Symbol(name, "_return"), WI = 0.0, simple_well = false, type = :btes, kwarg...)
   
    return supply_well, return_well

end

struct BTESWellSupplyToReturnMassCT <: Jutul.AdditiveCrossTerm
    btes_cell::Int64
end

struct BTESWellSupplyToReturnEnergyCT <: Jutul.AdditiveCrossTerm
    btes_cell::Int64
end