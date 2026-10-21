# Included into KAPreconditioners by isolated probes. Stable arguments at this
# boundary also avoid capturing reassigned per-level loop variables.
@inline function perf_level_interval(rows, first, count, worker, n, workers)
    low_row = 1 + (worker - 1)*n ÷ workers
    high_row = worker*n ÷ workers
    last = first + count - 1
    low = searchsortedfirst(rows, low_row, first, last, Base.Order.Forward)
    high = searchsortedlast(rows, high_row, first, last, Base.Order.Forward)
    return low, high
end

@noinline function perf_level_rows!(f::F, rows, first, count, n, minbatch, fixed) where F
    workers = Threads.nthreads(:default)
    if fixed && count >= 2*workers
        @batch for worker in 1:workers
            low, high = perf_level_interval(rows, first, count, worker, n, workers)
            @inbounds for q in low:high
                f(rows[q])
            end
        end
    elseif fixed
        @inbounds for q in first:(first + count - 1)
            f(rows[q])
        end
    else
        @batch minbatch = minbatch for q in 1:count
            @inbounds f(rows[first + q - 1])
        end
    end
    return nothing
end
