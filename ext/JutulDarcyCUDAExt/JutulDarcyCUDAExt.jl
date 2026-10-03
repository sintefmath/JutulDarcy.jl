module JutulDarcyCUDAExt
    using Jutul, JutulDarcy, CUDA, LinearAlgebra, SparseArrays, Adapt
    import Jutul: @tic

    # The host model owns well names, domains and range maps. CUDA kernels only
    # need the merged well count from this metadata.
    Adapt.adapt_structure(::CUDA.KernelAdaptor, info::JutulDarcy.MultiWellInfo) =
        JutulDarcy.KernelMultiWellInfo(length(info))

    # MultiSegmentWell's symbolic type and name are host metadata. Its
    # numerical fields remain available to device kernels.
    function Adapt.adapt_structure(to::CUDA.KernelAdaptor,
            well::JutulDarcy.MultiSegmentWell)
        return JutulDarcy.MultiSegmentWell(
            nothing, well.num_nodes, well.num_segments, well.num_perforations,
            Adapt.adapt(to, well.perforations),
            Adapt.adapt(to, well.neighborship),
            Adapt.adapt(to, well.end_nodes),
            Adapt.adapt(to, well.surface),
            nothing,
            Adapt.adapt(to, well.segment_models),
            Adapt.adapt(to, well.multiwell)
        )
    end

    JutulDarcy.kernel_abstractions_backend(::Val{:ka_cuda}) = CUDA.CUDABackend()

    timeit_debug_enabled() = Jutul.timeit_debug_enabled()

    include("ilu0.jl")
    include("krylov.jl")
    include("cuda_utils.jl")
    include("cpr.jl")
end
