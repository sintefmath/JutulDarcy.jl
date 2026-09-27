# GPU support

JutulDarcy includes support for running the simulator as kernels. This means that the entire simulator is massively parallel, which is highly beneficial for running cases on GPU. or on CPUs with a large number of threads. Everything can run on GPU including property evaluation, equation assembly and linear solves. This is by far the fastest way to run simulations that have more than about 100,000 cells.

## How to use

If you have installed JutulDarcy, you should start by adding the desired backend that mathces your GPU device:

### Adding CUDA backend

```julia
using Pkg
Pkg.add("CUDA") # Requires a CUDA-capable GPU
```

### Adding AMDGPU backend

```julia
using Pkg
Pkg.add("AMDGPU") # Requires a AMD ROC-capable GPU
```

### Loading a model

Once the packages have been added to the same environment as JutulDarcy, you can load them to enable GPU support. Let us grab the first ten steps of the EGG benchmark model:

```julia
using Jutul, JutulDarcy
dpth = JutulDarcy.GeoEnergyIO.test_input_file_path("EGG", "EGG.DATA")
case = setup_case_from_data_file(dpth)
case = case[1:10]
```

### Running on CPU

If we wanted to run this in parallel on CPU we would simply call `simulate_reservoir` with `:ka` as the mode:

```julia
result_cpu = simulate_reservoir(case, mode = :ka);
```

### Running on CUDA GPUs

```julia
using CUDA
result_cpu = simulate_reservoir(case, mode = :ka_cuda);
```

### Running on CUDA GPUs

```julia
using CUDA
result_cpu = simulate_reservoir(case, mode = :ka_cuda);
```

### Additional options

The default settings use a fairly conservative setup strategy for the AMG and the other parts of the linear solver. For instance, many cases can run more quickly if we reuse AMG operators and use the vendor ILU(0):

```julia
upd = :operators # Reuse operators, can be very fast, but degrades performance
simulate_reservoir(case;
    linear_solver_arg = (
        update_type = upd,
        smoother_type = :vendor_ilu
    )
)
```
