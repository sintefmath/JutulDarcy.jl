module JutulDarcyMetalExt
    using JutulDarcy, Metal

    JutulDarcy.kernel_abstractions_backend(::Val{:ka_metal}) = Metal.MetalBackend()
end
