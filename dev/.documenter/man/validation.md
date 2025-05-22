
# Overview {#Overview}

JutulDarcy.jl uses industry standard discretizations and will match other simulators when the parameters and driving forces are identical. Setting up identical parameters can be challenging, however, as there are many subtleties in exactly how input is interpreted. In this section, we go over some things to keep in mind when doing a direct comparison to other simulators.

## Comparison of simulations {#Comparison-of-simulations}

If you want to compare two numerical models it is advisable to start with a simple physical process that you fully understand how to set up in both simulators. In this section, we provide some hints to help this process.

### Input files {#Input-files}

The easiest way to perform a one-to-one comparison is through input files (.DATA files). JutulDarcy takes care to interpret these files in a similar manner to commercial codes, including features like transmissibility multipliers, block-to-block transmissibilities and processing of corner-point meshes around faults and eroded layers. Take note of yellow text during model setup - this can indicate unsupported or partially supported features. Some of these messages may refer to features that have little or no impact on simulation outcome, and others may have substantial impacts. For a direct comparison we recommend removing the keywords that give warnings and run the modified files in both JutulDarcy and the comparison simulator. Note that different simulators can sometimes also handle defaulted keywords differently, even when both simulators come from the same vendor. If you are experiencing mismatch, make sure that the case is not dependent on any defaulted values.

JutulDarcy does in general not make use of keywords from input files that describe the solution strategy like linear or nonlinear tolerance adjustments, tuning of time-stepping, switching to IMPES/AIM, etc.

### Time-steps {#Time-steps}

JutulDarcy performs dynamic timestepping, taking care to hit the provided report steps to provide output. For large report steps, different simulators can take different time-steps which results in differences in results. Some simulators, like Eclipse, also report sparse data like well results at a higher interval than the requested reporting interval which can make direct comparison of plots tricky. To get higher resolution in the same manner for JutulDarcy, setting `output_substates = true` is useful.

### Tolerances {#Tolerances}

The default tolerances in JutulDarcy are fairly strict, with measures for both point-wise and integral errors. Direct comparison may require checking that convergence criteria are set to similar values.

### Parallelism {#Parallelism}

If runtimes are to be compared, make sure that the simulators are using a similar model for parallelization. The fastest solves are generally achieved with either MPI or GPU parallelization, both of which require a bit of extra effort to set up.

## Downloadable examples {#Downloadable-examples}

The remainder of this part of the manual are a set of examples where JutulDarcy is compared against other simulators that are widely in use (e.g. OPM Flow, MRST, AD-GPRS, Eclipse 100/300). Extensive validation has also been performed on non-public models that for obvious reason cannot be directly incorporated in the documentation. All these examples automatically download the prerequisite data together with one or more results to compare with. The examples should be directly reproducible, but may require additional installation/licenses to produce the comparison results.

### User examples {#User-examples}

If you have an example that can be shared, either for the documentation or for regression testing, please open an issue. If you have a confidential model with a known solution that can be shared privately, feel free to [reach out directly](mailto:olav.moyner@sintef.no).
