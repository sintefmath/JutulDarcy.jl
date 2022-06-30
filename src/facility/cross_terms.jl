struct ReservoirWellCrossTerm{T<:AbstractVector, I<:AbstractVector} <: Jutul.AdditiveCrossTerm
    WI::T
    reservoir_cells::I
    well_cells::I
end
