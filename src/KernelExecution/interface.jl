"""
    kernel_abstractions_backend(mode)

Return the KernelAbstractions backend associated with a reservoir simulator
mode. The KernelExecution module implements the CPU modes and package
extensions add accelerator-specific modes.
"""
function kernel_abstractions_backend end

kernel_abstractions_backend(mode::Symbol) =
    kernel_abstractions_backend(Val(mode))

function kernel_abstractions_backend(::Val{mode}) where mode
    throw(ArgumentError(
        "No KernelAbstractions backend is available for mode :$mode. " *
        "Load the corresponding backend package or pass ka_backend explicitly."))
end
