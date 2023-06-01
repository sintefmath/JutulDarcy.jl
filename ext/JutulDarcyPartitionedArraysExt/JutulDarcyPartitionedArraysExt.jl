module JutulDarcyPartitionedArraysExt
    using Jutul, JutulDarcy
    # Specific dependencies
    using PartitionedArrays, MPI
    # Already in Jutul
    # using SparseArrays, Krylov, LinearAlgebra#, LinearOperators

    import Jutul: PArraySimulator, MPISimulator, PArrayExecutor
    import Jutul: DebugPArrayBackend, JuliaPArrayBackend, MPI_PArrayBackend
    import Jutul: partition_distributed, simulate_parray
    import JutulDarcy: reservoir_partition, partitioner_input

    function JutulDarcy.simulate_reservoir_parray(case::JutulCase, np::Int, backend::PArrayBackend; conn = :unit)
        N, T, groups = partitioner_input(case.model, case.parameters, conn = conn)
        rmodel = reservoir_model(case.model)
        nc = number_of_cells(rmodel.domain)
        p_num = partition_distributed(N, T, nc = nc, np = np, groups = groups)
        p = reservoir_partition(case.model, p_num)
        return simulate_parray(case, p, backend)
    end

end
