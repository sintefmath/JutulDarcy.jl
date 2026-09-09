module JutulDarcyAMDGPUExt
    using JutulDarcy, AMDGPU

    function __init__()
        JutulDarcy.register_kernel_abstractions_backend!(
            Symbol("ka-amd"), () -> AMDGPU.ROCBackend())
    end
end
