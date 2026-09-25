module JutulDarcyAMDGPUExt
    using JutulDarcy, AMDGPU

    JutulDarcy.kernel_abstractions_backend(::Val{:ka_amd}) = AMDGPU.ROCBackend()
end
