using Jutul, JutulDarcy, MultiComponentFlash, Test

@testset "GASWAT deck physics" begin
    deck = """
RUNSPEC
METRIC
GASWAT
WATER
DIFFUSE
COMPS
3 /
DIMENS
3 1 1 /
GRID
DXV
3*100 /
DYV
100 /
DZV
20 /
TOPS
3*1100 /
PORO
3*0.2 /
PERMX
3*100 /
PERMY
3*100 /
PERMZ
3*10 /
PROPS
EOS
PR /
CNAMES
H2O C1 CO2 /
MW
18.015 16.043 44.01 /
TCRIT
647.3 190.6 304.7 /
PCRIT
220.483 46.0421 73.8659 /
VCRIT
0.056 0.098 0.094 /
ACF
0.344 0.013 0.225 /
BIC
0.072 0.003732278 0.01 /
DENSITY
1* 1020.3 /
PVTW
215 1.0223 4.1483E-5 0.32929 0 /
WSF
0.2 0
0.95 1 /
GSF
0 0 0.1
0.05 0 0.1
0.79 0.95 0.8
0.8 1 0.8 /
TEMPVD
1210 60 /
ZMFVD
1210 0.004511 0.970489 0.025 /
DIFFCGAS
0.026 0.013 0.052 /
SOLUTION
EQUIL
1210 1* 1210 0.1 2* 2* 0 2 /
SCHEDULE
FBHPDEF
50 1000 /
WELSPECS
WELL FIELD 1 1 1210 GAS /
/
COMPDAT
WELL 1 1 1 1 OPEN 2* 0.2 /
/
WELLSTRE
GAS 0 0.8 0.2 /
/
WINJGAS
WELL STREAM GAS /
/
WCONINJE
WELL GAS OPEN BHP 2* 90 /
/
TSTEP
0.1 /
WELLSHUT
WELL /
/
TSTEP
0.1 /
WCONPROD
WELL OPEN GRAT 2* 1000 /
/
WELOPEN
WELL OPEN /
/
TSTEP
0.1 /
"""
    mktempdir() do dir
        filename = joinpath(dir, "GASWAT.DATA")
        write(filename, deck)
        data = JutulDarcy.GeoEnergyIO.parse_data_file(filename)
        @test data["RUNSPEC"]["DIFFUSE"]
        @test data["PROPS"]["DIFFCGAS"] ≈ [0.026, 0.013, 0.052]./86400
        @test data["PROPS"]["GSF"][1][end, 3] ≈ 8e4
        raw = JutulDarcy.GeoEnergyIO.parse_data_file(filename; units = nothing)
        @test raw["SCHEDULE"]["STEPS"][1]["FBHPDEF"] == [50.0, 1000.0]
        case = setup_case_from_data_file(data; extra_outputs = true)
        model = reservoir_model(case)
        sys = model.system
        @test number_of_phases(sys) == 2
        @test JutulDarcy.number_of_components(sys) == 3
        @test sys.equation_of_state.type isa SoreideWhitson
        @test !haskey(model.primary_variables, :ImmiscibleSaturation)
        @test haskey(case.parameters[:Reservoir], :VaporDiffusivities)
        @test haskey(model.data_domain, :vapor_diffusion, Cells())
        @test case.forces[end][:Facility].limits[:WELL].bhp ≈ 5e6
        # Diffusion acts directly on flash mole fractions. Reverse the face
        # orientation and remove gas to check conservation and phase cutoff.
        mw = MultiComponentFlash.molar_masses(sys.equation_of_state)
        x = [1.0, 0.0, 0.0]
        yl, yr = [0.01, 0.89, 0.10], [0.01, 0.79, 0.20]
        flash(y) = FlashedMixture2Phase(MultiComponentFlash.two_phase_lv,
            ones(3), 0.5, x, y, 1.0, 1.0)
        diffusion_state = (
            FlashResults = [flash(yl), flash(yr)],
            VaporDiffusivities = case.parameters[:Reservoir][:VaporDiffusivities],
            PhaseMassDensities = [1000.0 1000.0; 5.0 7.0],
            Saturations = [0.2 0.4; 0.8 0.6])
        q0 = JutulDarcy.SVector{3}(0.0, 0.0, 0.0)
        function flux(gradient, state = diffusion_state)
            Dl = JutulDarcy.phase_diffusivities(state, LiquidPhase())
            Dv = JutulDarcy.phase_diffusivities(state, VaporPhase())
            return JutulDarcy.add_diffusive_component_flux(q0, Dl, Dv,
                1, state, model, gradient, Val(3))
        end
        q = flux(TPFA(1, 2, 1))
        @test abs(q[1]) < 1e-20
        @test q[2] > 0.0 && q[3] < 0.0
        concentration = (5.0/sum(mw.*yl) + 7.0/sum(mw.*yr))/2
        expected = -diffusion_state.VaporDiffusivities[:, 1].*0.6.*concentration.*mw.*(yr - yl)
        @test q ≈ expected
        @test flux(TPFA(2, 1, -1)) ≈ -q
        @test all(iszero, flux(TPFA(1, 2, 1),
            merge(diffusion_state, (Saturations = [0.2 1.0; 0.8 0.0],))))
        initial = case.state0[:Reservoir]
        @test all(isfinite, initial[:Pressure])
        @test maximum(abs, sum(initial[:OverallMoleFractions]; dims = 1) .- 1) < 1e-10
        # Connate water must survive conversion from hydrostatic phase volumes
        # to the overall composition used by the compositional flash.
        @test initial[:Saturations][1, :] ≈ fill(0.2, 3) atol = 1e-6
        for c in axes(initial[:OverallMoleFractions], 2)
            p = initial[:Pressure][c]
            T = case.parameters[:Reservoir][:Temperature][c]
            z = initial[:OverallMoleFractions][:, c]
            f = MultiComponentFlash.flashed_mixture_2ph(sys.equation_of_state,
                (p = p, T = T, z = z))
            sl, sv = MultiComponentFlash.phase_saturations(sys.equation_of_state, p, T, f)
            @test [sl, sv] ≈ initial[:Saturations][:, c] atol = 1e-6
            @test f.vapor.mole_fractions[2]/f.vapor.mole_fractions[3] ≈ 0.970489/0.025 rtol = 1e-7
        end
        result = simulate_reservoir(case; info_level = -1,
            initial_dt = 864.0, tol_cnv = 1e-5)
        @test length(result.states) == 3
        @test all(s -> all(isfinite, s[:Pressure]), result.states)
        @test all(s -> maximum(abs, sum(s[:Saturations]; dims = 1) .- 1) < 1e-8, result.states)
        # Phase-prefixed component arrays must also reach the kernel flux.
        kernel_result = simulate_reservoir(case; mode = :ka, linear_solver = nothing,
            info_level = -1, initial_dt = 864.0, tol_cnv = 1e-5)
        @test length(kernel_result.states) == length(result.states)
        for (state, reference) in zip(kernel_result.states, result.states)
            @test state[:Pressure] ≈ reference[:Pressure] rtol = 1e-7
            @test state[:OverallMoleFractions] ≈ reference[:OverallMoleFractions] rtol = 1e-7
        end
    end
end

@testset "Søreide–Whitson flash at vapor endpoint" begin
    # A hydrogen-rich production-well mixture can pass the stability test but
    # converge to V = 1 in SSI. It must have single-phase properties and AD.
    mw = [18.015, 16.043, 28.013, 44.01, 30.07, 44.097, 58.124, 2.016]./1000
    pc = [220.483, 46.0421, 33.9439, 73.8659, 48.8387, 42.4552, 37.9665, 13].*1e5
    tc = [647.3, 190.6, 126.2, 304.7, 305.43, 369.8, 425.2, 33.2]
    vc = [0.056, 0.098, 0.09, 0.094, 0.148, 0.2, 0.255, 0.065]./1000
    acf = [0.344, 0.013, 0.04, 0.225, 0.0986, 0.1524, 0.201, -0.218]
    bic = [0.072, 0, 0.00373227833333338, 0.01, 0, 0.01, 0.01, 0, 0.01,
        0, 0.01, 0, 0.01, zeros(9)..., -0.00748975999999986, 0, 0.01, 0.01, 0.01, 0]
    aij = zeros(8, 8)
    k = 1
    for i in 2:8, j in 1:i-1
        aij[i, j] = aij[j, i] = bic[k]
        k += 1
    end
    mixture = MultiComponentMixture(MolecularProperty.(mw, pc, tc, vc, acf);
        names = ["H2O", "C1", "N2", "CO2", "C2", "C3", "NC4", "H2"], A_ij = aij)
    eos = GenericCubicEOS(mixture, SoreideWhitson(mixture);
        volume_shift = (0.22, 0, 0, 0, 0, -0.08349872, 0, -0.2565496))
    options = JutulDarcy.FlashResults{SSIFlash, false, false}(
        SSIFlash(), false, 1e-8, 10.0, false, false)
    initial = JutulDarcy.static_flashed_mixture(eos, Float64)
    u = [9.424, 0.0026, 0.1747, 0.01813, 0.005029, 0.0014356, 7.973e-5, 3.985e-5]
    composition(u) = JutulDarcy.SVector{8}(u[2:8]..., 1-sum(u[2:8]))
    cond = (p = 1e6*u[1], T = 333.15, z = composition(u))
    v, _, _ = JutulDarcy.numeric_flash(initial, options, eos, cond)
    @test v == 1.0
    flash(u) = JutulDarcy.immutable_flash_result(initial, options, eos,
        1e6*u[1], 333.15, composition(u), 0.0)
    result = flash(u)
    @test result.state == MultiComponentFlash.single_phase_v
    @test result.V == 1.0
    @test result.vapor.mole_fractions == composition(u)
    function properties(u)
        result = flash(u)
        sat = MultiComponentFlash.phase_saturations(eos, 1e6*u[1], 333.15, result)
        rho = MultiComponentFlash.mass_densities(eos, 1e6*u[1], 333.15, result)
        return [sat.S_l, rho[1], rho[2]]
    end
    ad = JutulDarcy.ForwardDiff.jacobian(properties, u)
    fd = similar(ad)
    for i in eachindex(u)
        if i == 1
            h = 1e-5
        else
            h = 1e-8
        end
        plus, minus = copy(u), copy(u)
        plus[i] += h
        minus[i] -= h
        fd[:, i] = (properties(plus)-properties(minus))/(2h)
    end
    @test all(iszero, ad[1, :])
    @test ad ≈ fd rtol = 1e-6
end
