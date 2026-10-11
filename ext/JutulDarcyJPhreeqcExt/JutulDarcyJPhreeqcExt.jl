module JutulDarcyJPhreeqcExt
    using Jutul, JutulDarcy, JPhreeqc
    import ForwardDiff
    import JutulDarcy.Tracers as Tracers
    import JutulDarcy: setup_phreeqc_model!, phreeqc_component_names,
        phreeqc_initial_concentrations, phreeqc_component_molar_masses, close_phreeqc!

    include("types.jl")
    include("chemistry.jl")
    include("variables.jl")
    include("api.jl")
end
