module JutulDarcyMetalExt
    using JutulDarcy, Metal

    function __init__()
        JutulDarcy.register_kernel_abstractions_backend!(
            Symbol("ka-metal"), () -> Metal.MetalBackend())
    end
end
