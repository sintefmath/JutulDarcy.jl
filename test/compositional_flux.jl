using Jutul, JutulDarcy, StaticArrays, Test

@testset "Compositional flux sparsity tracing" begin
    liquid = [0.2 0.4; 0.3 0.1; 0.5 0.5]
    vapor = [0.6 0.1; 0.1 0.7; 0.3 0.2]
    trace_fractions(X) = Jutul.SparsityTracingWrapper(
        allocate_array_ad(X; npartials = 3, tag = Cells()))
    liquid_traced = trace_fractions(liquid)
    vapor_traced = trace_fractions(vapor)
    upw = SPU(1, 2)

    @testset "Advection, aqueous=$aqueous" for aqueous in (false, true)
        n = aqueous ? 4 : 3
        phase_indices = aqueous ? (1, 2, 3) : (1, 2)
        mass_fluxes = aqueous ? (7.0, 2.0, -3.0) : (2.0, -3.0)
        traced_fluxes = aqueous ? (
            Jutul.as_tracer(7.0, 5),
            Jutul.as_tracer(2.0, 3),
            Jutul.as_tracer(-3.0, 4)
        ) : (Jutul.as_tracer(2.0, 3), Jutul.as_tracer(-3.0, 4))
        q = zeros(SVector{n, Float64})
        q_traced = zeros(MVector{n, typeof(first(traced_fluxes))})
        numeric = JutulDarcy.compositional_fluxes!(q, liquid, vapor,
            mass_fluxes, upw, Val(aqueous), phase_indices, Val(3))
        traced = JutulDarcy.compositional_fluxes!(q_traced,
            liquid_traced, vapor_traced, traced_fluxes, upw,
            Val(aqueous), phase_indices, Val(3))

        expected = 2.0*liquid[:, 1] - 3.0*vapor[:, 2]
        aqueous && push!(expected, 7.0)
        @test numeric ≈ expected
        @test traced isa typeof(q_traced)
        @test map(x -> x.val, traced) ≈ expected
        for i in 1:3
            @test Jutul.ST.deriv(traced[i]).nzind == [1, 2, 3, 4]
        end
        if aqueous
            @test Jutul.ST.deriv(traced[4]).nzind == [5]
        end
    end

    @testset "Diffusion" begin
        masses = [0.018, 0.016, 0.002]
        D = reshape([1e-9, 2e-9, 4e-9], 3, 1)
        saturations = [0.2 0.4; 0.8 0.6]
        densities = [1000.0 1100.0; 5.0 7.0]
        grad = TPFA(1, 2, 1)
        q = zeros(SVector{4, Float64})
        q_traced = zeros(MVector{4, typeof(liquid_traced[1, 1])})
        # The fourth component belongs to the aqueous phase and must be kept.
        q = setindex(q, 7.0, 4)
        q_traced = setindex(q_traced, Jutul.as_tracer(7.0, 3), 4)
        numeric = JutulDarcy.add_phase_diffusive_component_flux(q, D, 1,
            grad, saturations, densities, liquid, 1, masses, Val(3))
        traced = JutulDarcy.add_phase_diffusive_component_flux(q_traced, D,
            1, grad, saturations, densities, liquid_traced, 1, masses, Val(3))

        @test traced isa typeof(q_traced)
        @test map(x -> x.val, traced) ≈ numeric
        for i in 1:3
            @test Jutul.ST.deriv(traced[i]).nzind == [1, 2]
        end
        @test traced[4].val == 7.0
        @test Jutul.ST.deriv(traced[4]).nzind == [3]
    end
end
