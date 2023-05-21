module JutulDarcyHYPREExt
    using Jutul, JutulDarcy, HYPRE

    function JutulDarcy.update_pressure_system!(A_p::HYPRE.HYPREMatrix, A, w_p, bz, ctx)
        error()
        T_p = eltype(A_p)
        nz = nonzeros(A_p)
        nz_s = nonzeros(A)
        cols = Jutul.colvals(A)
        @assert size(nz) == size(nz_s)
        n = size(A_p, 1)
        # Update the pressure system with the same pattern in-place
        tb = minbatch(ctx, n)
        for row in 1:n
            update_row_csr!(nz, A_p, w_p, cols, nz_s, row)
        end
    end
end
