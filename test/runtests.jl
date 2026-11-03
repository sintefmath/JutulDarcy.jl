using JutulDarcy, TestItems, Test, TestItemRunner

@testitem "Mesh and reservoir" tags = [:ci_core, :mesh, :input, :wells] begin
    include("mesh_and_reservoir.jl")
    include("corner_geometry.jl")
    include("schedule_completions.jl")
    include("zero_rate_deck_controls.jl")
end

@testitem "Utilities" tags = [:ci_core, :utilities, :input, :wells] begin
    include("utils.jl")
    include("secondary_scheme.jl")
    include("input_timestep_defaults.jl")
    include("efficiency_summary.jl")
    include("group_summary.jl")
end

@testitem "Rel. Perm." tags = [:ci_core, :properties, :relperm, :initialization] begin
    include("relperm.jl")
    include("initialization_scaling.jl")
end

@testitem "Capillary pressure" tags = [:ci_core, :properties, :capillary] begin
    include("capillary.jl")
end

@testitem "Single-phase" tags = [:ci_flow, :flow, :singlephase] begin
    include("singlephase.jl")
end

@testitem "KernelAbstractions single model" tags = [:ci_backends, :backends, :kernel_abstractions, :gpu, :flow] begin
    include("kernel_abstractions.jl")
end

@testitem "CPRW BHP controls" tags = [:ci_solvers, :solvers, :cpr, :wells] begin
    include("cprw_bhp.jl")
end

@testitem "Float32 reservoir models" tags = [:ci_backends, :backends, :float32, :flow] begin
    include("reservoir_float32.jl")
end

@testitem "Multi-phase" tags = [:ci_flow, :flow, :immiscible, :sources] begin
    include("multiphase.jl")
    include("source_terms.jl")
end

@testitem "Multi-model (wells + reservoir)" tags = [:ci_wells, :flow, :wells, :blackoil, :compositional, :immiscible] begin
    include("multimodel.jl")
end

@testitem "Closed-loop cross terms (deprecated)" tags = [:ci_wells, :wells, :closed_loop] begin
    include("closed_loop_crossterms.jl")
end

@testitem "Runtime well shutting" tags = [:ci_wells, :wells, :controls] begin
    include("well_shutting.jl")
end

@testitem "Well group controls" tags = [:ci_wells, :wells, :controls, :groups] begin
    include("well_groups.jl")
end

@testitem "Well limits" tags = [:ci_wells, :wells, :controls] begin
    include("well_limits.jl")
end

@testitem "Multisegment wells" tags = [:ci_wells, :wells, :multisegment] begin
    include("mswells.jl")
end

@testitem "Sensitivities (simple)" tags = [:ci_adjoint_parameters, :adjoint, :parameters, :immiscible] begin
    include("sens_bl.jl")
end

@testitem "MRST input cases" tags = [:ci_flow, :input, :blackoil, :compositional, :wells] begin
    include("mrst_cases.jl")
end

@testitem "PArray solve" tags = [:ci_partitioned, :parallel, :mpi, :solvers] begin
    # JLD2 can resolve OrderedCollections to TestItemServer's vendored copy.
    # Keep disk round trips in a clean process with the same environment,
    # bounds checks and coverage flags as this test worker.
    script = joinpath(@__DIR__, "parray.jl")
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())`
    if Base.JLOptions().code_coverage != 0
        coverage_file = joinpath(@__DIR__, "..", "parray-lcov.info")
        cmd = `$cmd --code-coverage=$coverage_file`
    end
    @test success(pipeline(`$cmd $script`; stdout = stdout, stderr = stderr))
end

@testitem "Scalarization" tags = [:ci_core, :backends, :scalarization] begin
    include("scalarization.jl")
end

@testitem "Multi-phase systems basics" tags = [:ci_core, :properties, :systems] begin
    include("systems.jl")
end

@testitem "Thermal" tags = [:ci_thermal, :thermal, :flow, :wells] begin
    include("thermal.jl")
end

@testitem "Thermal DATA and passive tracers" tags = [:ci_compositional, :thermal, :input, :tracers] begin
    include("thermal_data.jl")
end

@testitem "Tabulated thermal properties" tags = [:ci_thermal, :thermal, :properties, :tabulation] begin
    include("thermal_properties.jl")
end

@testitem "Tabulated thermal backends" tags = [:ci_compositional, :thermal, :properties, :tabulation, :backends, :kernel_abstractions, :gpu] begin
    include("thermal_backends.jl")
end

@testitem "NLDD" tags = [:ci_solvers, :solvers, :nldd, :parallel] begin
    include("nldd.jl")
end

@testitem "Sensitivities (multimodel)" tags = [:ci_adjoint_parameters, :adjoint, :parameters, :optimization, :wells, :blackoil, :compositional, :immiscible] begin
    include("sens_multimodel.jl")
end

@testitem "Sensitivities (forces)" tags = [:ci_adjoint_controls, :adjoint, :controls, :wells] begin
    include("sens_forces.jl")
end

@testitem "Sensitivities (state0)" tags = [:ci_adjoint_parameters, :adjoint, :initialization] begin
    include("sens_state0.jl")
end

@testitem "Sensitivities (equil)" tags = [:ci_adjoint_parameters, :adjoint, :initialization, :equilibration] begin
    include("sens_equil.jl")
end

@testitem "Compositional validation" tags = [:ci_compositional, :compositional, :flash, :flow] begin
    include("compositional.jl")
end

@testitem "Compositional flux sparsity tracing" tags = [:ci_compositional, :compositional, :flux, :sparsity] begin
    include("compositional_flux.jl")
end

@testitem "Component diffusion parameters and fluxes" tags = [:ci_compositional, :compositional, :diffusion, :flux] begin
    include("diffusion.jl")
end

@testitem "GASWAT deck physics" tags = [:ci_compositional, :compositional, :co2, :input] begin
    include("gas_water.jl")
end

@testitem "Two-phase flash derivatives" tags = [:ci_compositional, :compositional, :flash, :derivatives] begin
    include("flash_implicit_derivatives.jl")
end

@testitem "CO2 props" tags = [:ci_compositional, :co2, :properties] begin
    include("co2props.jl")
end

@testitem "Discretizations" tags = [:ci_flow, :discretization, :input, :blackoil] begin
    include("discretizations.jl")
end

@testitem "Sequential schemes" tags = [:ci_core, :solvers, :sequential] begin
    include("sequential.jl")
end

@testitem "GPU" tags = [:ci_compositional, :gpu, :backends, :compositional, :flash] begin
    include("gpu.jl")
end

@testitem "Well constructors" tags = [:ci_wells, :wells, :mesh] begin
    include("well_constructor.jl")
end

@testitem "Merged wells" tags = [:ci_core, :wells, :mesh] begin
    include("merge_wells.jl")
end

@testitem "HistoryMatching" tags = [:ci_core, :history_matching, :optimization] begin
    include("historymatch.jl")
end

@testitem "DFM" tags = [:ci_solvers, :fractures, :dfm, :discretization] begin
    include("dfm.jl")
end

@testitem "Sensitivities (inner adjoint)" tags = [:ci_adjoint_controls, :adjoint, :parameters, :optimization, :wells, :blackoil, :compositional, :immiscible] begin
    include("sens_inner_adjoint.jl")
end

@testitem "CI group coverage" tags = [:ci_core, :infrastructure] begin
    using TestItemRunner
    include("ci_groups.jl")
    checked = Ref(0)
    TestItemRunner.run_tests(joinpath(@__DIR__, ".."); filter = ti -> begin
        testitem_in_ci_group(ti)
        checked[] += 1
        return false
    end)
    @test checked[] > 0
end

include("ci_groups.jl")
@run_package_tests filter = testitem_in_ci_group verbose = true
nothing
