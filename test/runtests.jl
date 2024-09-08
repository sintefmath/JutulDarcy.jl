using JutulDarcy, TestItems, Test, TestItemRunner

@testitem "Utilities" begin
    include("utils.jl")
end

@testitem "Rel. Perm." begin
    include("relperm.jl")
end

@testitem "Single-phase" begin
    include("singlephase.jl")
end

@testitem "Multi-phase" begin
    include("multiphase.jl")
end

@testitem "Multi-model (wells + reservoir)" begin
    include("multimodel.jl")
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

@testitem "NLDD" begin
    include("nldd.jl")
end

@testitem "Sensitivities (multimodel)" begin
    include("sens_multimodel.jl")
end

@testitem "Compositional validation" begin
    include("compositional.jl")
end

@run_package_tests
nothing
# include("gpu.jl")
