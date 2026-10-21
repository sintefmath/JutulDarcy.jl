# Optional Linux benchmark affinity, separate from the production library.
function configure_benchmark_pinning!()
    requested = get(ENV, "PERF_PIN_CPUS", "")
    isempty(requested) && return
    @assert Sys.islinux() "PERF_PIN_CPUS requires Linux"
    cpus = parse.(Int, split(requested, ','))
    @assert length(cpus) == Threads.nthreads(:default)
    observed = similar(cpus)
    Threads.@threads :static for slot in eachindex(cpus)
        cpu = cpus[slot]
        @assert 0 <= cpu < 1024
        mask = zeros(UInt8, 128)
        mask[1 + cpu ÷ 8] = UInt8(1) << (cpu % 8)
        status = ccall(:sched_setaffinity, Cint, (Cint, Csize_t, Ptr{UInt8}), 0, length(mask), mask)
        @assert status == 0 "Could not set worker affinity"
        observed[slot] = ccall(:sched_getcpu, Cint, ())
    end
    @assert observed == cpus
    println("Pinned default workers to CPUs ", observed)
    flush(stdout)
end
configure_benchmark_pinning!()
