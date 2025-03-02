using Test, JutulDarcy

import JutulDarcy.CO2Properties: compute_co2_brine_props

@testset "CO2-Brine Properties" begin
    p = 250e5
    T = 273.0 + 30

    props = compute_co2_brine_props(p, T, [0.0], ["NaCl"])
    @test props == compute_co2_brine_props(p, T)

    @test props[:K] ≈ [0.00447321, 37.0489] rtol = 1e-3
    @test props[:viscosity] ≈ [0.000800893, 0.000106415] rtol = 1e-3
    @test props[:density] ≈ [1006.68, 913.948] rtol = 1e-3

    p = 180e5
    T = 273.0 + 60
    props = compute_co2_brine_props(180e5, T)

    @test props[:K] ≈ [0.00843427, 46.5601] rtol = 1e-3
    @test props[:viscosity] ≈ [0.000470367, 5.96663e-5] rtol = 1e-3
    @test props[:density] ≈ [991.11, 668.644] rtol = 1e-3

    props = compute_co2_brine_props(180e5, T, [0.05], ["NaCl"])
    @test props[:K] ≈ [0.00843427, 88.0887] rtol = 1e-3
    @test props[:viscosity] ≈ [0.00064919, 5.96663e-5] rtol = 1e-3
    @test props[:density] ≈ [1092.47, 668.644] rtol = 1e-3

    props = compute_co2_brine_props(180e5, T, [0.01, 0.01, 0.005, 0.01, 0.01, 0.012], ["NaCl", "KCl", "CaSO4", "CaCl2", "MgSO4", "MgCl2"])
    @test props[:K] ≈ [0.00843427, 143.368] rtol = 1e-3
    @test props[:viscosity] ≈ [0.000501, 5.9666e-5] rtol = 1e-3
    @test props[:density] ≈ [1170.37, 668.644] rtol = 1e-3
end
