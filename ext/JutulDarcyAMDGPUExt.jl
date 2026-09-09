module JutulDarcyAMDGPUExt
    using JutulDarcy, AMDGPU

    function __init__()
        JutulDarcy.register_kernel_abstractions_backend!(
            :ka_amd, () -> AMDGPU.ROCBackend())
    end
end
