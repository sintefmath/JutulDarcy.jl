using JutulDarcy, TestItems, Test, TestItemRunner

@testitem "Mesh and reservoir" begin
    include("mesh_and_reservoir.jl")
    include("corner_geometry.jl")
    include("schedule_completions.jl")
    include("zero_rate_deck_controls.jl")
end

@testitem "Utilities" begin
    include("utils.jl")
    include("secondary_scheme.jl")
    include("input_timestep_defaults.jl")
    include("efficiency_summary.jl")
end

@testitem "Rel. Perm." begin
    include("relperm.jl")
    include("initialization_scaling.jl")
end

@testitem "Capillary pressure" begin
    include("capillary.jl")
end

@testitem "Single-phase" begin
    include("singlephase.jl")
end

@testitem "KernelAbstractions single model" begin
    include("kernel_abstractions.jl")
end

@testitem "CPRW BHP controls" begin
    include("cprw_bhp.jl")
end

@testitem "Float32 reservoir models" begin
    include("reservoir_float32.jl")
end

@testitem "Multi-phase" begin
    include("multiphase.jl")
    include("source_terms.jl")
end

@testitem "Multi-model (wells + reservoir)" begin
    include("multimodel.jl")
end

@testitem "Runtime well shutting" begin
    include("well_shutting.jl")
end

@testitem "Multisegment wells" begin
    include("mswells.jl")
end

@testitem "Sensitivities (simple)" begin
    include("sens_bl.jl")
end

@testitem "MRST input cases" begin
    include("mrst_cases.jl")
end

@testitem "PArray solve" begin
    include("parray.jl")
end

@testitem "Scalarization" begin
    include("scalarization.jl")
end

@testitem "Multi-phase systems basics" begin
    include("systems.jl")
end

@testitem "Thermal" begin
    include("thermal.jl")
end

@testitem "Tabulated thermal properties" begin
    include("thermal_properties.jl")
end

@testitem "Tabulated thermal backends" begin
    include("thermal_backends.jl")
end

@testitem "NLDD" begin
    include("nldd.jl")
end

@testitem "Sensitivities (multimodel)" begin
    include("sens_multimodel.jl")
end

@testitem "Sensitivities (forces)" begin
    include("sens_forces.jl")
end

@testitem "Sensitivities (state0)" begin
    include("sens_state0.jl")
end

@testitem "Sensitivities (equil)" begin
    include("sens_equil.jl")
end

@testitem "Compositional validation" begin
    include("compositional.jl")
end

@testitem "Compositional flux sparsity tracing" begin
    include("compositional_flux.jl")
end

@testitem "Component diffusion parameters and fluxes" begin
    include("diffusion.jl")
end

@testitem "GASWAT deck physics" begin
    include("gas_water.jl")
end

@testitem "Two-phase flash derivatives" begin
    include("flash_implicit_derivatives.jl")
    include("flash_reuse.jl")
end

@testitem "CO2 props" begin
    include("co2props.jl")
end

@testitem "Discretizations" begin
    include("discretizations.jl")
end

@testitem "Sequential schemes" begin
    include("sequential.jl")
end

@testitem "GPU" begin
    include("gpu.jl")
end

@testitem "Well constructors" begin
    include("well_constructor.jl")
end

@testitem "Merged wells" begin
    include("merge_wells.jl")
end

@testitem "HistoryMatching" begin
    include("historymatch.jl")
end

@testitem "DFM" begin
    include("dfm.jl")
end

@run_package_tests
nothing
