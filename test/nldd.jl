using JutulDarcy, Test, HYPRE, MPI, PartitionedArrays

@testset "NLDD/MPI" begin
    spe1_pth = JutulDarcy.GeoEnergyIO.test_input_file_path("spe1", "BENCH_SPE1.DATA")
    il = -1
    case = setup_case_from_data_file(f, split_wells = true);
    ws_base, = simulate_reservoir(case, method = :newton, info_level = il);

    ws = Dict()
    ws[:nldd], = simulate_reservoir(case, method = :nldd, info_level = il);
    ws[:nldd_4], = simulate_reservoir(case, method = :nldd, info_level = il, nldd_arg = Dict(:no_blocks => 4));

    ws[:aspen], = simulate_reservoir(case, method = :aspen, info_level = il);
    ws[:mpi], = simulate_reservoir(case, method = :newton, mode = :mpi, info_level = il);
    ws[:nldd_mpi], = simulate_reservoir(case, method = :nldd, mode = :mpi, info_level = il);

    get_prod_rate(x) = x[:PRODUCER][Symbol("Surface total rate")]
    get_grat(x) = x[:PRODUCER][Symbol("Surface gas rate")]
    get_prod_bhp(x) = x[:PRODUCER][Symbol("Bottom hole pressure")]

    for (k, ws_compare) in pairs(ws)
        @testset "$k producer" begin
            rate_ref = get_prod_rate(ws_base)
            rate = get_prod_rate(ws_compare)
            @test rate_ref ≈ rate rtol = 1e-3

            rate_ref = get_grat(ws_base)
            rate = get_grat(ws_compare)
            @test rate_ref ≈ rate rtol = 1e-3

            pbhp_ref = get_prod_bhp(ws_base)
            pbhp = get_prod_bhp(ws_compare)
            @test pbhp_ref ≈ pbhp rtol = 1e-2
        end
    end
end