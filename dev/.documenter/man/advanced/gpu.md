
# GPU support {#GPU-support}

JutulDarcy includes experimental support for running linear solves on the GPU. For many simulations, the linear systems are the most compute-intensive part and a natural choice for acceleration. At the moment, the support is limited to CUDA GPUs through [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl). For the most efficient CPR preconditioner, [AMGX.jl](https://github.com/JuliaGPU/AMGX.jl) is required which is currently limited to Linux systems. Windows users may have luck by running Julia inside [WSL](https://learn.microsoft.com/en-us/windows/wsl/install).

## How to use {#How-to-use}

If you have installed JutulDarcy, you should start by adding the CUDA and optionally the AMGX packages using the package manager:

```julia
using Pkg
Pkg.add("CUDA") # Requires a CUDA-capable GPU
Pkg.add("AMGX") # Requires CUDA + Linux
```


Once the packages have been added to the same environment as JutulDarcy, you can load them to enable GPU support. Let us grab the first ten steps of the EGG benchmark model:

```julia
using Jutul, JutulDarcy
dpth = JutulDarcy.GeoEnergyIO.test_input_file_path("EGG", "EGG.DATA")
case = setup_case_from_data_file(dpth)
case = case[1:10]
```


### Running on CPU {#Running-on-CPU}

If we wanted to run this on CPU we would simply call `simulate_reservoir`:

```julia
result_cpu = simulate_reservoir(case);
```


### Running on GPU with block ILU(0) {#Running-on-GPU-with-block-ILU0}

If we now load `CUDA` we can run the same simulation using the CUDA-accelerated linear solver. By itself, CUDA only supports the ILU(0) preconditioner. JutulDarcy will automatically pick this preconditioner when CUDA is requested without AMGX, but we write it explicitly here:

```julia
using CUDA
result_ilu0_cuda = simulate_reservoir(case, linear_solver_backend = :cuda, precond = :ilu0);
```


### Running on GPU with CPR AMGX-ILU(0) {#Running-on-GPU-with-CPR-AMGX-ILU0}

Loading the AMGX package makes a pure GPU-based two-stage CPR available. Again, we are explicit in requesting CPR, but if both `CUDA` and `AMGX` are available and functional this is redundant:

```julia
using AMGX
result_amgx_cuda = simulate_reservoir(case, linear_solver_backend = :cuda, precond = :cpr);
```


In short, load `AMGX` and `CUDA` and run `simulate_reservoir(case, linear_solver_backend = :cuda)` to get GPU results. The EGG model is quite small, so if you want to see significant performance increases, a larger case will be necessary. `AMGX` also contains a large number of options that can be configured for advanced users.

## Technical details and limitations {#Technical-details-and-limitations}

The GPU implementation relies on assembly on CPU and pinned memory to transfer onto the CPU. This means that the performance can be significantly improved by launching Julia with multiple threads to speed up the non-GPU parts of the code. AMGX is currently single-GPU only and does not work with MPI. To make use of lower precision, specify `Float32` in the `float_type` argument to the linear solver. Additional arguments to `AMGX` can also be specified this way. For example, we can solve using aggregation AMG in single precision by doing the following:

```julia
simulate_reservoir(case,
    linear_solver_backend = :cuda,
    linear_solver_arg = (
        float_type = Float32,
        algorithm = "AGGREGATION",
        selector = "SIZE_8"
        )
    )
```


::: warning Experimental status

Multiple successive runs with different `AMGX` instances have resulted in crashes when old instances are garbage collected. This part of the code is still considered experimental, with contributions welcome if you are using it.

:::
