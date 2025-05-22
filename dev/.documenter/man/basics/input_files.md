
# Input formats {#Input-formats}

It is also possible to read cases that have been set up in MRST (see [`setup_case_from_mrst`](/man/basics/input_files#JutulDarcy.setup_case_from_mrst) and [`simulate_mrst_case`](/man/basics/input_files#JutulDarcy.simulate_mrst_case)) or from .DATA files (see [`setup_case_from_data_file`](/man/basics/input_files#JutulDarcy.setup_case_from_data_file) and [`simulate_data_file`](/man/basics/input_files#JutulDarcy.simulate_data_file))

## MAT-files from the Matlab Reservoir Simulation Toolbox (MRST) {#MAT-files-from-the-Matlab-Reservoir-Simulation-Toolbox-MRST}

### Simulation of .MAT files {#Simulation-of-.MAT-files}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.setup_case_from_mrst' href='#JutulDarcy.setup_case_from_mrst'><span class="jlbinding">JutulDarcy.setup_case_from_mrst</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
setup_case_from_mrst("filename.mat"; kwarg...)
```


Set up a [`Jutul.JutulCase`](/ref/jutul#Jutul.JutulCase) from a MRST-exported .mat file.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/input_simulation/mrst_input.jl#L1102-L1106" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.simulate_mrst_case' href='#JutulDarcy.simulate_mrst_case'><span class="jlbinding">JutulDarcy.simulate_mrst_case</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
ws, states = simulate_mrst_case(file_name)
simulate_mrst_case(file_name; <keyword arguments>)
```


Simulate a MRST case from `file_name` as exported by `writeJutulInput` in MRST.

**Arguments**
- `file_name::String`: The path to a `.mat` or `.data` file that is to be simulated.
  

**Keyword arguments**
- `extra_outputs::Vector{Symbol} = [:Saturations]`: Additional variables to output from the simulation.
  
- `write_output::Bool = true`: Write output (in the default JLD2 format)
  
- `output_path = nothing`: Directory for output files. Files will be written under this directory. Defaults to the folder of `file_name`.
  
- `write_mrst = true`: Write MRST compatible output after completed simulation that can be read by `readJutulOutput` in MRST.
  
- `backend=:csc`: choice of backend for linear systems. `:csc` for default Julia sparse, `:csr` for experimental parallel CSR.
  
- `verbose=true`: print some extra information specific to this routine upon calling
  
- `nthreads=Threads.nthreads()`: number of threads to use
  
- `linear_solver=:bicgstab`: name of Krylov.jl solver to use, or :direct (for small cases only)
  
- `info_level=0`: standard Jutul info_level. 0 for minimal printing, -1 for no printing, 1-5 for various levels of verbosity
  

Additional input arguments are passed onto, [`setup_case_from_mrst`](/man/basics/input_files#JutulDarcy.setup_case_from_mrst), [`setup_reservoir_simulator`](/man/highlevel#JutulDarcy.setup_reservoir_simulator) and [`simulator_config`](/ref/jutul#Jutul.simulator_config-Tuple{Any}) if applicable.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/input_simulation/mrst_input.jl#L1494-L1525" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### MRST-specific types {#MRST-specific-types}
<details class='jldocstring custom-block' open>
<summary><a id='Jutul.MRSTWrapMesh' href='#Jutul.MRSTWrapMesh'><span class="jlbinding">Jutul.MRSTWrapMesh</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MRSTWrapMesh(G, N = nothing)
```


Mesh that adapts an exported MRST mesh to the Jutul interface. `G` is assumed to be read directly from file using `MAT.matread`. The raw exported grid can be found under the `data` field.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/Jutul.jl/blob/v0.3.11/src/meshes/mrst.jl#L9-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## DATA-files from commercial reservoir modelling software {#DATA-files-from-commercial-reservoir-modelling-software}

JutulDarcy can set up cases from Eclipse-type input files by making use of the [GeoEnergyIO.jl](https://github.com/sintefmath/GeoEnergyIO.jl) package for parsing. This package is a direct dependency of JutulDarcy and these cases can be simulated directly. If you want to parse the input files and possibly modify them in your Julia session before the case is simulated, we refer you to the [GeoEnergyIO Documentation](https://sintefmath.github.io/GeoEnergyIO.jl/dev/).

If you want to directly simulate a file from disk, you can sue the high level functions that automatically parse the files for you:
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.simulate_data_file' href='#JutulDarcy.simulate_data_file'><span class="jlbinding">JutulDarcy.simulate_data_file</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
simulate_data_file(inp; parse_arg = NamedTuple(), kwarg...)
```


Simulate standard input file (with extension .DATA). `inp` can either be the output from the GeoEnergyIO function `parse_data_file` or a `String` for the path of an input file with the .DATA extension.

Additional arguments are passed onto [`simulate_reservoir`](/man/highlevel#JutulDarcy.simulate_reservoir). Extra inputs to the parser can be sent as a `setup_arg` `NamedTuple`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/input_simulation/data_input.jl#L1-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.setup_case_from_data_file' href='#JutulDarcy.setup_case_from_data_file'><span class="jlbinding">JutulDarcy.setup_case_from_data_file</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
case = setup_case_from_data_file(
    filename;
    parse_arg = NamedTuple(),
    kwarg...
)

case, data = setup_case_from_data_file(filename; include_data = true)
```


Set up a [`JutulCase`](/ref/jutul#Jutul.JutulCase) from a standard input file (with extension .DATA). Additional arguments to the parser can be passed as key/values in a `NamedTuple` given as `parse_arg`. The optional input `include_data=true` will make the function return the parsed data in addition the case itself. Additional keyword arguments are passed on to [`setup_case_from_parsed_data`](/man/basics/input_files#JutulDarcy.setup_case_from_parsed_data).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/input_simulation/data_input.jl#L21-L35" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.setup_case_from_parsed_data' href='#JutulDarcy.setup_case_from_parsed_data'><span class="jlbinding">JutulDarcy.setup_case_from_parsed_data</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
setup_case_from_parsed_data(datafile; skip_wells = false, simple_well = true, use_ijk_trans = true, verbose = false, kwarg...)
```


Set up a case from a parsed input file (in `Dict` format). Internal function, not exported. Use [`setup_case_from_data_file`](/man/basics/input_files#JutulDarcy.setup_case_from_data_file).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/input_simulation/data_input.jl#L57-L62" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.convert_co2store_to_co2_brine' href='#JutulDarcy.convert_co2store_to_co2_brine'><span class="jlbinding">JutulDarcy.convert_co2store_to_co2_brine</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
convert_co2store_to_co2_brine(data; verbose = true)
```


Converts a CO2STORE data file to a co2-brine model. The conversion should be close to equivialent for models without salt. The data will be copied if modifications are necessary.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/input_simulation/converter.jl#L1-L7" target="_blank" rel="noreferrer">source</a></Badge>

</details>


Reservoir simulator input files are highly complex and contain a great number of different keywords. JutulDarcy and GeoEnergyIO has extensive support for this format, but many keywords are missing or only partially supported. Sometimes cases can be made to run by removing keywords that do not impact simulation outcomes, or have very little impact. If you want a turnkey open source solution for simulation reservoir models you should also have a look at [OPM Flow](https://opm-project.org/).

If you can share your input file or the missing keywords [in the issues section](https://github.com/sintefmath/JutulDarcy.jl/issues) it may be easier to figure out if a case can be supported.
