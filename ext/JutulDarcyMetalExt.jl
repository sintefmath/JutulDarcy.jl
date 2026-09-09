module JutulDarcyMetalExt
    using JutulDarcy, Metal

    function __init__()
        JutulDarcy.register_kernel_abstractions_backend!(
            :ka_metal, () -> Metal.MetalBackend())
    end
end
