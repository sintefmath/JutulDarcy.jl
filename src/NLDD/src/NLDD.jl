module NLDD
    export NLDDSimulator, partition_uniform_1d
    import Jutul:   update_before_step!,
                    update_after_step!,
                    solve_ministep, perform_step!,
                    store_output!,
                    reset_state_to_previous_state!,
                    simulator_config,
                    final_simulation_message

    using Jutul, Base.Threads
    using JutulDarcy
    using Polyester
    using DataStructures
    using LinearAlgebra
    using TimerOutputs
    using ProgressMeter
    using ForwardDiff
    using SparseArrays, StaticArrays
    using Printf
    import Jutul: @tic

    # All types in module
    include("types.jl")
    # DD simulator
    include("simulator.jl")
    include("transfer.jl")
    include("utils.jl")
    include("reservoir.jl")
    # Adaptivity
    include("adaptive.jl")
    include("logger.jl")
    # ASPEN
    include("partial_update.jl")
    include("aspen.jl")
    include("aspen_assembly.jl")
    # Tests etc
    include("test_functions.jl")
end
