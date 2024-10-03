function coarsen_reservoir(D::DataDomain, partition; functions = Dict())
    if !haskey(functions, :permeability)
        functions[:permeability] = Jutul.CoarsenByHarmonicAverage()
    end
    if !haskey(functions, :porosity)
        functions[:porosity] = Jutul.CoarsenByVolumeAverage()
    end
    return Jutul.coarsen_data_domain(D, partition, functions = functions)
end
