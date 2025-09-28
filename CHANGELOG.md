# JutulDarcy.jl Changelog

This changelog documents all notable changes to JutulDarcy.jl based on the GitHub releases.

## [v0.2.46] - 2024-07-06

### New Features

- Better treatment of initial state for black-oil models in sensitivities code
- Improved examples
- Improved performance of flux calculations
- Several improvements to multisegment well setup (allowing for example variable roughness along wellbore)
- Multisegment wells should converge more easily in the very low flow rate regime (for example for a disabled well)
- Geothermal models now use tables at constant pressure instead of volume for heat capacity
- Well QoIs in objective function now support symbol aliases - can specify e.g. `:lrat` instead of providing a type
- Additional well setup routines are now flexible with respect to numeric types, making AD usage easier
- Updated minimum Jutul version: Latest version improves accuracy of adjoint gradients for models with many substeps and optimization problems with significant drift in model parameters during optimization

### Merged Pull Requests

- Put sequential code in module + updates ([#130](https://github.com/sintefmath/JutulDarcy.jl/pull/130))
- Type fix to BlackOil variable ([#137](https://github.com/sintefmath/JutulDarcy.jl/pull/137))
- Improve gradients example, improve initial state gradients ([#138](https://github.com/sintefmath/JutulDarcy.jl/pull/138))
- Multisegment well updates ([#142](https://github.com/sintefmath/JutulDarcy.jl/pull/142))
- Fix face average inbounds call ([#143](https://github.com/sintefmath/JutulDarcy.jl/pull/143))
- Correct phase properties and indexing in multicomponent bc ([#146](https://github.com/sintefmath/JutulDarcy.jl/pull/146))
- Update differentiable geothermal example ([#149](https://github.com/sintefmath/JutulDarcy.jl/pull/149))
- MSwell improvements + geothermal updates ([#150](https://github.com/sintefmath/JutulDarcy.jl/pull/150))
- Fix objective format in compositional example ([#151](https://github.com/sintefmath/JutulDarcy.jl/pull/151))

### Closed Issues

- Minor bug in multicomponent/bc.jl ([#145](https://github.com/sintefmath/JutulDarcy.jl/pull/145))

## [v0.2.45] - 2024-06-17

## [v0.2.44] - 2024-04-09

### New Features

- New example: [Advanced history matching](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/data_assimilation/advanced_history_match)
- Fix to tracers for non-standard discretizations (e.g. WENO/AvgMPFA)
- Additional options for net-present-value objective
- Brooks-Corey parametric relative permeability now has configurable fields
- Improvements to alternative timesteppers on control changes

### Merged Pull Requests

- Create DocPreviewCleanup.yml ([#113](https://github.com/sintefmath/JutulDarcy.jl/pull/113))
- Timestepping updates ([#114](https://github.com/sintefmath/JutulDarcy.jl/pull/114))
- Add example demonstrating blending to automatically discover regions ([#115](https://github.com/sintefmath/JutulDarcy.jl/pull/115))
- Updates to gradients ([#116](https://github.com/sintefmath/JutulDarcy.jl/pull/116))

### Closed Issues

- turn off print ([#109](https://github.com/sintefmath/JutulDarcy.jl/pull/109))

## [v0.2.43] - 2024-03-19

### New Features

- Support for gradients with respect to well controls/limits, forces (see [new example for NPV optimization](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/workflow/rate_optimization))
- Fix to radius for friction model in MS wells when using convenience constructors
- New timestepper that reduces timestep on large changes in wells (preview, not exported)

### Merged Pull Requests

- Add minimum timestep to simulator config ([#99](https://github.com/sintefmath/JutulDarcy.jl/pull/99))
- Support for gradients of well controls/composition/temperature, boundary conditions and source terms ([#105](https://github.com/sintefmath/JutulDarcy.jl/pull/105))
- Timestepping ([#106](https://github.com/sintefmath/JutulDarcy.jl/pull/106))
- Make sure segment radius is set to wellbore radius multisegment wells ([#107](https://github.com/sintefmath/JutulDarcy.jl/pull/107))
- Rate optimization utilities + example ([#108](https://github.com/sintefmath/JutulDarcy.jl/pull/108))

### Closed Issues

- StackOverflowError ([#104](https://github.com/sintefmath/JutulDarcy.jl/pull/104))

## [v0.2.42] - 2024-03-07

### New Features

- Salinity is now accounted for in brine viscosity for CO2-brine model
- New ease-of-use feature for hydrostatic equilibrium (example: [https://sintefmath.github.io/JutulDarcy.jl/dev/examples/workflow/equilibrium_state](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/workflow/equilibrium_state))
- Bugfixes to setup of thermal models
- Bugfixes to initialization of equilibrium pressure for aquifer models
- Highly experimental support for writing SUMMARY type output
- Temperature based thresholding in NLDD
- Equilibriation now supports more density functions
- Fix to irreversible rock compressibility

### Merged Pull Requests

- Increase swap in CI ([#98](https://github.com/sintefmath/JutulDarcy.jl/pull/98))
- Add utility and example for making hydrostatic equilibrium initialization easier in scripts ([#101](https://github.com/sintefmath/JutulDarcy.jl/pull/101))
- Add salinity to brine viscosity in CO2 model, fix to irreversible rock compressibility, bugfixes to aquifer initialization, faster setup of models ([#102](https://github.com/sintefmath/JutulDarcy.jl/pull/102))
- Make sure rock and component heat capacities are set to DataDomain ([#103](https://github.com/sintefmath/JutulDarcy.jl/pull/103))

### Closed Issues

- objective function modification ([#100](https://github.com/sintefmath/JutulDarcy.jl/pull/100))

## [v0.2.41] - 2024-02-14

### Major New Features

- **Expanded geothermal features**, including borehole thermal energy storage (BTES) support and automatic setup of property tables with support for salinity
- **Support for single phase and multiphase tracers** that can have per-tracer custom behavior and interact with other properties of the system
- **A validated polymer model** (viscosity changes, dead pore space, permeability reduction, adsorption) that makes use of the new tracer functionality
- **Routines to rebuild an existing reservoir model** with for example new wells or changed properties. This is useful when you have loaded a complex model from input files and want to alter it without changing the input file
- **Docs have been restructured and revised**

### New Examples

#### Workflow examples

- [Adding tracers to a flow simulation](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/workflow/tracers_two_wells)
- [Adding new wells to an existing model](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/workflow/adding_new_wells)

#### Validation examples

- [Polymer injection in a 2D black-oil reservoir model](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/validation/validation_polymer)
- [Aquifer thermal energy storage (ATES) validation](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/validation/validation_thermal)

#### Geothermal examples

- [Borehole Thermal Energy Storage (BTES)](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/geothermal/btes)
- [Geothermal doublet](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/geothermal/geothermal_doublet)
- [High-temperature Aquifer Thermal Energy Storage (HT-ATES)](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/geothermal/htates_intro)

### Improvements

- Improved performance for some BLAS versions
- Examples have been reorganized by topic, both in the code and in the documentation to make navigation easier
- Many fixes and improvements to thermal well setup, including additional parameters exposed for sensitivity analysis
- Support for SLGOF keyword
- Bumped Makie compatibility

### Merged Pull Requests

- Support for different types of tracers + simple polymer model ([#90](https://github.com/sintefmath/JutulDarcy.jl/pull/90))
- Reorganize the examples a bit ([#91](https://github.com/sintefmath/JutulDarcy.jl/pull/91))
- Improve polymer model + tracers ([#92](https://github.com/sintefmath/JutulDarcy.jl/pull/92))
- Updated docs, added support for tracers + polymer and improved equilibriation ([#93](https://github.com/sintefmath/JutulDarcy.jl/pull/93))
- First steps towards geothermal module ([#94](https://github.com/sintefmath/JutulDarcy.jl/pull/94))
- Geothermal improvements ([#96](https://github.com/sintefmath/JutulDarcy.jl/pull/96))
- Enable BTES example ([#97](https://github.com/sintefmath/JutulDarcy.jl/pull/97))

## [v0.2.40] - 2024-01-16

### Merged Pull Requests

- Update tests, robustness for init ([#89](https://github.com/sintefmath/JutulDarcy.jl/pull/89))

## [v0.2.39] - 2024-01-08

### Breaking Changes

- **Wells are no longer multisegment wells by default**
- **Changes in the inner `compute_well_qoi` interface** (now takes in well state and facility state instead of global state)

### Other Updates

- Function for getting Brooks-Corey capillary pressure (`brooks_corey_pc`)
- Doc improvements
- Improvements to accuracy and robustness for compositional models
- Support for VCRIT as alternative to ZCRIT
- Support for storing phase fluxes (`store_total_fluxes`/`store_phase_fluxes`). Not exported, name subject to change.

### Merged Pull Requests

- Update doc blocks in example ([#82](https://github.com/sintefmath/JutulDarcy.jl/pull/82))
- Doc updates 3 ([#83](https://github.com/sintefmath/JutulDarcy.jl/pull/83))
- Convert warnings to errors in docs, update docs and simplify compute_well_qoi, ZCRIT support ([#84](https://github.com/sintefmath/JutulDarcy.jl/pull/84))
- Compositional improvements ([#85](https://github.com/sintefmath/JutulDarcy.jl/pull/85))

## [v0.2.38] - 2024-11-28

### New Features

- **Improved closed loop geothermal support** (including reinjection support)
- **New differentiability interface** that makes it easy to optimize w.r.t. any parameter used in model setup
- **Improved plotting behavior on MacOS** for complex 3D plots

The new functionality is demonstrated in two new examples:
- [A fully differentiable geothermal doublet: History matching and control optimization](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/workflow/fully_differentiable_geothermal)
- [Adjoint gradients for the SPE1 model](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/data_assimilation/spe1_gradients)

In addition, the [Buckley-Leverett adjoint example](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/data_assimilation/optimize_simple_bl) has been revised to include this new functionality.

### Merged Pull Requests

- Fix typos ([#121](https://github.com/sintefmath/JutulDarcy.jl/pull/121))
- Make test deterministic and less flaky ([#122](https://github.com/sintefmath/JutulDarcy.jl/pull/122))
- Update vitepress deploy ([#123](https://github.com/sintefmath/JutulDarcy.jl/pull/123))
- Geothermal improvements ([#124](https://github.com/sintefmath/JutulDarcy.jl/pull/124))
- Support for Jutul 0.4 + type fixes for generic differentiation ([#125](https://github.com/sintefmath/JutulDarcy.jl/pull/125))
- Add test for state0 sensitivities ([#128](https://github.com/sintefmath/JutulDarcy.jl/pull/128))
- Update adjoint for well controls to latest Jutul + fixes ([#129](https://github.com/sintefmath/JutulDarcy.jl/pull/129))
- New optimization interface+ example ([#136](https://github.com/sintefmath/JutulDarcy.jl/pull/136))

### Closed Issues

- Using Clapeyron.jl for properties? ([#135](https://github.com/sintefmath/JutulDarcy.jl/pull/135))

## [v0.2.37] - 2024-11-16

### New Features

- **Support for WENO schemes in all solvers**
- **Two new examples** demonstrating consistent schemes (NTPFA, AvgMPFA) and high-resolution (WENO)
- **CPR support for alternative discretizations**
- **Refactored linear solver internal interface**
- **Improved boundary transmissibility calculations** for flow bc
- **Support for timestepping based on composition change**

### Merged Pull Requests

- CPRW support ([#78](https://github.com/sintefmath/JutulDarcy.jl/pull/78))
- AMGX finalizer fix + update to Jutul 0.30 ([#79](https://github.com/sintefmath/JutulDarcy.jl/pull/79))
- WENO support, improved CPR support for non-standard discretizations, more tests ([#81](https://github.com/sintefmath/JutulDarcy.jl/pull/81))

## [v0.2.36] - 2024-11-13

### New Features

- **Support for linear solves (CPR/ILU) on GPU for CUDA**. [See docs for more details.](https://sintefmath.github.io/JutulDarcy.jl/dev/man/advanced/gpu)
- **Improved EDIT and TRANX/TRANY/TRANZ support**
- **Bug fixes in CPR**
- **Significantly improved performance for adaptive NLDD solvers**
- **Example updates**

### Merged Pull Requests

- Adding example of hybrid simulation with DNN for rel perm ([#75](https://github.com/sintefmath/JutulDarcy.jl/pull/75))
- Tentative CUDA/AMGX support ([#76](https://github.com/sintefmath/JutulDarcy.jl/pull/76))
- Add high-level CUDA support, doc improvements, NLDD performance increases, CSR backend as default ([#77](https://github.com/sintefmath/JutulDarcy.jl/pull/77))

### Closed Issues

- setting custom objective function ([#69](https://github.com/sintefmath/JutulDarcy.jl/pull/69))
- DisableControl() and compute_well_qoi ([#73](https://github.com/sintefmath/JutulDarcy.jl/pull/73))
- `extra_out` error for `setup_reservoir_model()` ([#74](https://github.com/sintefmath/JutulDarcy.jl/pull/74))

## [v0.2.35] - 2024-10-15

### New Features

- **Improved docs**
- **Coarsening/upscaling routines**
- **Parametric rel.perm. curves (LET/Brooks-Corey)**
- **Updated plotting, new plotter for field-scale responses**
- **Units and lots of improvements in well plotter**
- **Support for SALINITY keyword in CO2 models** as an alternative to using salt mole fractions

### New Examples

- Added a new example demonstrating history matching a coarse model using the CGNet method. This includes setting up the optimization problem and comparing results with a fine model.
- Introduced a new example showing how to coarsen a model using various methods and comparing the results with the fine-scale model.
- New example demonstrating different relative permeability functions and how these can be calibrated against data

### Merged Pull Requests

- Added `step_no` to `sum_co2_inventory` ([#68](https://github.com/sintefmath/JutulDarcy.jl/pull/68))
- Upscaling/coarsening of models, CGNet example, relperm example, LET functions, new GUIs, rewritten well plotter ([#71](https://github.com/sintefmath/JutulDarcy.jl/pull/71))

### Closed Issues

- bottom hole pressure ([#63](https://github.com/sintefmath/JutulDarcy.jl/pull/63))
- Function arguments in examples ([#67](https://github.com/sintefmath/JutulDarcy.jl/pull/67))
- setting custom objective function ([#69](https://github.com/sintefmath/JutulDarcy.jl/pull/69))

## [v0.2.34] - 2024-10-01

### Merged Pull Requests

- Rock hysteresis, fix scaled capillary pressure regression, bump HYPRE ([#66](https://github.com/sintefmath/JutulDarcy.jl/pull/66))

## [v0.2.33] - 2024-09-28

### New Features

- **Improved docs**
- **Revised PINCH support**
- **A limited version of initial state gradients**
- **New example of CO2 property calculations**

### Merged Pull Requests

- Documentation updates ([#65](https://github.com/sintefmath/JutulDarcy.jl/pull/65))

### Closed Issues

- multiple well gradients ([#60](https://github.com/sintefmath/JutulDarcy.jl/pull/60))

## [v0.2.32] - 2024-09-24

### New Features

- **Salinity support for CO2-brine model** (port of code from [@lsalo](https://github.com/lsalo))
- **Support for older Julia versions** (@gbruer15)
- **New docs**

### Merged Pull Requests

- Switch docs over to Vitepress ([#59](https://github.com/sintefmath/JutulDarcy.jl/pull/59))
- Updated DelimitedFiles compat ([#61](https://github.com/sintefmath/JutulDarcy.jl/pull/61)) (@gbruer15)
- Support for salts in CO2 model, new property codes for CO2 and brine based on work by Lluis Salo ([#62](https://github.com/sintefmath/JutulDarcy.jl/pull/62))

### Closed Issues

- Support for Julia 1.7 and 1.8 ([#50](https://github.com/sintefmath/JutulDarcy.jl/pull/50))

## [v0.2.31] - 2024-09-23

### Updates

- **Support for numerical aquifers**
- **Better validation of static reservoir properties**
- **Cleaned up capillary pressure definitions for two-phase scenarios**
- **Fix `generate_jutuldarcy_examples` read-only issue**
- **Improved docs**

### Merged Pull Requests

- Support for numerical aquifers and cleanup of capillary pressure ([#58](https://github.com/sintefmath/JutulDarcy.jl/pull/58))

## [v0.2.30] - 2024-09-08

### Merged Pull Requests

- Bugfix to initialization and fault plots ([#56](https://github.com/sintefmath/JutulDarcy.jl/pull/56))

## [v0.2.29] - 2024-09-08

### Merged Pull Requests

- New examples and example improvements ([#55](https://github.com/sintefmath/JutulDarcy.jl/pull/55))

## [v0.2.28] - 2024-08-31

### New Features

- **Simplified thermal implementation** (no more composite system usage)
- **Support for new consistent discretizations** with performant assembly: AvgMPFA and nonlinear TPFA. To use, pass `kgrad = :avgmpfa` or `kgrad = :ntpfa` to `setup_reservoir_model`.
- **Converter for CO2STORE models to K-value model**
- **MB convergence check support in compositional**
- **Fixes to initialization of water-gas models and thermal models**
- **Support for constant head aquifers**
- **Support for calculating well cells from trajectories**
- **Rewrote the `co2_sloped.jl` example** to demonstrate hysteresis, well trajectories and boundary conditions

### Merged Pull Requests

- Nonlinear finite-volume + AvgMPFA support ([#52](https://github.com/sintefmath/JutulDarcy.jl/pull/52))
- Simplify thermal implementation ([#53](https://github.com/sintefmath/JutulDarcy.jl/pull/53))
- Nonlinear finite-volume schemes, constant pressure aquifers, converter for CO2STORE, simplify thermal ([#54](https://github.com/sintefmath/JutulDarcy.jl/pull/54))

## [v0.2.27] - 2024-08-02

### New Features

- **Rewritten rel. perm. code**, should support all phase combinations
- **Hysteresis support for relative permeabilities** (Killough, Carlsen, Jargon and pure imbibition with drainage equilibriation). Compatible with endpoint scaling.
- **Support for injection at reservoir conditions (RESV-type)**
- **Improved docs & tests**

### Merged Pull Requests

- CompatHelper: bump compat for TestItems to 1, (keep existing compat) ([#48](https://github.com/sintefmath/JutulDarcy.jl/pull/48)) (@github-actions[bot])
- Hysteresis & relperm refactoring ([#51](https://github.com/sintefmath/JutulDarcy.jl/pull/51))

## [v0.2.26] - 2024-06-27

### Updates

- **Fixed crash in compositional initialization**
- **Fixed crash for two-point rel. perm. scaling for water-oil systems**
- **Greatly improved performance for immiscible viscosity evaluation**

## [v0.2.25] - 2024-06-25

### Updates

- **Support for WELTARG, WEFAC**
- **Improved initialization of complex models**
- **Better defaults for linear solver**

### Merged Pull Requests

- Improved performance, better behavior on input cases ([#47](https://github.com/sintefmath/JutulDarcy.jl/pull/47))

## [v0.2.24] - 2024-06-05

### Updates

- **Improvements to equilibriation** (maximum water saturation, non-unique data points for pressure, ...)
- **Better plotting of reservoirs and wells**
- **More accurate well indices** (effects from NTG + drainage radius)
- **Precompilation workflow now includes grid processing**
- **Support for AMGCL in CPR**
- **Better nonlinear performance for blackoil models**
- **Added examples for Eclipse-type input**: Norne, SPE1, SPE9, Egg and Olympus models with validation
- **Improved handling of many types of multipliers**
- **New unexported routine `reservoir_sensitivities`** computes sensitivities with respect to static properties (perm/poro/geometry/multipliers)
- **Output of thermal well results**
- **Better NLDD defaults**

### Merged Pull Requests

- Improved initialization, thermal output and grid processing ([#46](https://github.com/sintefmath/JutulDarcy.jl/pull/46))

### Closed Issues

- Error in simulate_mrst_case for firstJutulExample.m ([#38](https://github.com/sintefmath/JutulDarcy.jl/pull/38))
- Error in well output for JutulSPE1.m ([#39](https://github.com/sintefmath/JutulDarcy.jl/pull/39))

## [v0.2.23] - 2024-05-06

### Merged Pull Requests

- Improved adjoint support, additional grid processing keywords ([#42](https://github.com/sintefmath/JutulDarcy.jl/pull/42))

### Closed Issues

- Features Request ([#41](https://github.com/sintefmath/JutulDarcy.jl/pull/41))

## [v0.2.22] - 2024-04-14

### Changes

- **Allow `simulate_mrst_case` to support `.DATA` file format** (and write MRST output from that)
- **Fix issues running MRST exported cases for new MRST format** (issue [#39](https://github.com/sintefmath/JutulDarcy.jl/pull/39))

### Merged Pull Requests

- Fix MRST output writing ([#40](https://github.com/sintefmath/JutulDarcy.jl/pull/40))

## [v0.2.21] - 2024-04-12

### Changes

- **Documentation improvements**
- **Add support for stability bypass and flash reuse in compositional solvers** (pass `fast_flash=true` to `setup_reservoir_model` to enable)
- **Support for parsed volume shift in solvers**
- **Much faster CPR setup for many components** (no longer hard coded limit of up to 6 components for fastest version)
- **Improved compositional equilibriation support**
- **Well solver now is only active for first iteration when enabled by default**

### Merged Pull Requests

- Flash bypass and improved performance for many components ([#37](https://github.com/sintefmath/JutulDarcy.jl/pull/37))

## [v0.2.20] - 2024-04-04

### Merged Pull Requests

- Attempt at making GLMakie work in examples ([#34](https://github.com/sintefmath/JutulDarcy.jl/pull/34))
- Pre-solver for wells, improved MS well parsing, better NLDD support ([#35](https://github.com/sintefmath/JutulDarcy.jl/pull/35))
- Updates to docs + add CO2-brine model from CSP11 ([#36](https://github.com/sintefmath/JutulDarcy.jl/pull/36))

## [v0.2.19] - 2024-03-10

### Merged Pull Requests

- Add nonlinear domain decomposition solvers ([#33](https://github.com/sintefmath/JutulDarcy.jl/pull/33))

## [v0.2.18] - 2024-03-04

### Features

- **Support for parsing and setting up multisegment wells from .DATA input files**
- **Faster test suite**
- **Fix bug in capillary initialization when using EQUIL** - could lead to wrong sign for capillary pressure in subsequent simulations

#### Diffusion support in compositional models

- `simulate_reservoir` can now take a custom config object
- Option to specify well enthalpy during injection, either as function of pressure and temperature or a constant value
- New properties for brine-CO2 mixtures

#### Fixes

- Better handling of `ImmiscibleSystem` when a single phase is present
- K-value flash can now use analytical solution for two and three components
- Added scaled convergence criterion for thermal
- Single-phase equilibriation is now supported
- Heat capacity is now per component instead of per phase to support miscible systems

#### Keyword and input files

- Improve thermal support to include input files (WATDENT and keywords like WATVISCT)
- Handle WELOPEN keywords

### Merged Pull Requests

- Thermal rework ([#32](https://github.com/sintefmath/JutulDarcy.jl/pull/32))

## [v0.2.17] - 2024-02-07

### Changes

- **Support for PVT regions**
- **Support for WEFAC (well downtime)**
- **Remove parser and instead use GeoEnergyIO.jl package** where functionality has been moved
- **Plotting improvements**
- **Better initialization of various well keywords**
- **Initialization / equilibriation has been improved**

### Merged Pull Requests

- Remove parser that is moved to GeoEnergyIO ([#31](https://github.com/sintefmath/JutulDarcy.jl/pull/31))

## [v0.2.16] - 2024-01-28

### Merged Pull Requests

- K value flash support ([#30](https://github.com/sintefmath/JutulDarcy.jl/pull/30))

## [v0.2.15] - 2024-01-23

### New Features

- **Thermal support**
- **Improvements to plotting**
- **Parsing improvements (performance + keywords)**

## [v0.2.14] - 2024-01-08

### Merged Pull Requests

- Documentation updates ([#29](https://github.com/sintefmath/JutulDarcy.jl/pull/29))

## [v0.2.13] - 2024-01-04

### Merged Pull Requests

- Upgrades to .DATA and corner point grid parsers ([#28](https://github.com/sintefmath/JutulDarcy.jl/pull/28))

## [v0.2.12] - 2023-12-19

### New Features

- **Support for general corner point grids (GRDECL)**
- **Support for more keywords**
- **Serialization of primary variables over MPI**
- **Minor linear solver improvements in MPI mode**

### Merged Pull Requests

- Rewrite corner point support, improved parser, other small features ([#27](https://github.com/sintefmath/JutulDarcy.jl/pull/27))

## [v0.2.11] - 2023-11-17

### Merged Pull Requests

- CompatHelper: add new compat entry for Statistics at version 1, (keep existing compat) ([#23](https://github.com/sintefmath/JutulDarcy.jl/pull/23)) (@github-actions[bot])
- remove dependency of internal fields of `eos::AbstractEOS` ([#24](https://github.com/sintefmath/JutulDarcy.jl/pull/24))
- small clean-up of case setup from mrst data ([#25](https://github.com/sintefmath/JutulDarcy.jl/pull/25))
- CPR rewrite + performance improvements ([#26](https://github.com/sintefmath/JutulDarcy.jl/pull/26))

### Closed Issues

- CO2 brine example not working ([#22](https://github.com/sintefmath/JutulDarcy.jl/pull/22))

## [v0.2.10] - 2023-10-18

### Merged Pull Requests

- Better compositional solvers and fixes to simple wells + blackoil diffusion ([#21](https://github.com/sintefmath/JutulDarcy.jl/pull/21))

## [v0.2.9] - 2023-09-24

No specific changes documented.

## [v0.2.8] - 2023-08-07

### Merged Pull Requests

- Support for parsing and simulating .DATA files ([#19](https://github.com/sintefmath/JutulDarcy.jl/pull/19))

## [v0.2.7] - 2023-06-26

### Merged Pull Requests

- Port over JutulExamples cases ([#16](https://github.com/sintefmath/JutulDarcy.jl/pull/16))
- Add new ex to docs ([#17](https://github.com/sintefmath/JutulDarcy.jl/pull/17))
- Support for sensitivities for src/bc ([#18](https://github.com/sintefmath/JutulDarcy.jl/pull/18))

## [v0.2.6] - 2023-06-07

No specific changes documented.

## [v0.2.5] - 2023-06-06

No specific changes documented.

## [v0.2.4] - 2023-06-05

### Merged Pull Requests

- Migrate from SnoopPrecompile to PrecompileTools ([#14](https://github.com/sintefmath/JutulDarcy.jl/pull/14))
- MPI solve with HYPRE + PartitionedArrays ([#15](https://github.com/sintefmath/JutulDarcy.jl/pull/15))

## [v0.2.3] - 2023-04-22

### Merged Pull Requests

- Documentation updates and other small fixes ([#11](https://github.com/sintefmath/JutulDarcy.jl/pull/11))
- remove output path during precompilation ([#13](https://github.com/sintefmath/JutulDarcy.jl/pull/13))

### Closed Issues

- No writing access to `/tmp` leads to error during precompilation ([#12](https://github.com/sintefmath/JutulDarcy.jl/pull/12))

## [v0.2.2] - 2023-03-19

### Merged Pull Requests

- Add parameters kwarg to setup_reservoir_model ([#10](https://github.com/sintefmath/JutulDarcy.jl/pull/10)) (@gbruer15)

## [v0.2.1] - 2023-02-21

No specific changes documented.

## [v0.2.0] - 2023-02-04

Major version release with significant updates.

## [v0.1.6] - 2022-12-14

### New Features

- **Updates and compatibility fixes to MRST input**
- **Improved precompilation and tests**
- **RESV and historical RESV controls**
- **Can now run Norne field with a few missing features**

### Merged Pull Requests

- Support for VAPOIL ([#6](https://github.com/sintefmath/JutulDarcy.jl/pull/6))

### Merged Pull Requests

- Upgrades to .DATA and corner point grid parsers ([#28](https://github.com/sintefmath/JutulDarcy.jl/pull/28))

## [v0.2.12] - 2023-12-19

### New Features

- **Support for general corner point grids (GRDECL)**
- **Support for more keywords**
- **Serialization of primary variables over MPI**
- **Minor linear solver improvements in MPI mode**

### Merged Pull Requests

- Rewrite corner point support, improved parser, other small features ([#27](https://github.com/sintefmath/JutulDarcy.jl/pull/27))

## [v0.2.11] - 2023-11-17

### Merged Pull Requests

- CompatHelper: add new compat entry for Statistics at version 1, (keep existing compat) ([#23](https://github.com/sintefmath/JutulDarcy.jl/pull/23)) (@github-actions[bot])
- remove dependency of internal fields of `eos::AbstractEOS` ([#24](https://github.com/sintefmath/JutulDarcy.jl/pull/24))
- small clean-up of case setup from mrst data ([#25](https://github.com/sintefmath/JutulDarcy.jl/pull/25))
- CPR rewrite + performance improvements ([#26](https://github.com/sintefmath/JutulDarcy.jl/pull/26))

### Closed Issues

- CO2 brine example not working ([#22](https://github.com/sintefmath/JutulDarcy.jl/pull/22))

## [v0.2.10] - 2023-10-18

### Merged Pull Requests

- Better compositional solvers and fixes to simple wells + blackoil diffusion ([#21](https://github.com/sintefmath/JutulDarcy.jl/pull/21))

## [v0.2.9] - 2023-09-24

No specific changes documented.

## [v0.2.8] - 2023-08-07

### Merged Pull Requests

- Support for parsing and simulating .DATA files ([#19](https://github.com/sintefmath/JutulDarcy.jl/pull/19))

## [v0.2.7] - 2023-06-26

### Merged Pull Requests

- Port over JutulExamples cases ([#16](https://github.com/sintefmath/JutulDarcy.jl/pull/16))
- Add new ex to docs ([#17](https://github.com/sintefmath/JutulDarcy.jl/pull/17))
- Support for sensitivities for src/bc ([#18](https://github.com/sintefmath/JutulDarcy.jl/pull/18))

## [v0.2.6] - 2023-06-07

No specific changes documented.

## [v0.2.5] - 2023-06-06

No specific changes documented.

## [v0.2.4] - 2023-06-05

### Merged Pull Requests

- Migrate from SnoopPrecompile to PrecompileTools ([#14](https://github.com/sintefmath/JutulDarcy.jl/pull/14))
- MPI solve with HYPRE + PartitionedArrays ([#15](https://github.com/sintefmath/JutulDarcy.jl/pull/15))

## [v0.2.3] - 2023-04-22

### Merged Pull Requests

- Documentation updates and other small fixes ([#11](https://github.com/sintefmath/JutulDarcy.jl/pull/11))
- remove output path during precompilation ([#13](https://github.com/sintefmath/JutulDarcy.jl/pull/13))

### Closed Issues

- No writing access to `/tmp` leads to error during precompilation ([#12](https://github.com/sintefmath/JutulDarcy.jl/pull/12))

## [v0.2.2] - 2023-03-19

### Merged Pull Requests

- Add parameters kwarg to setup_reservoir_model ([#10](https://github.com/sintefmath/JutulDarcy.jl/pull/10)) (@gbruer15)

## [v0.2.1] - 2023-02-21

No specific changes documented.

## [v0.2.0] - 2023-02-04

Major version release with significant updates.

## [v0.1.6] - 2022-12-14

### New Features

- **Updates and compatibility fixes to MRST input**
- **Improved precompilation and tests**
- **RESV and historical RESV controls**
- **Can now run Norne field with a few missing features**

### Merged Pull Requests

- Support for VAPOIL ([#6](https://github.com/sintefmath/JutulDarcy.jl/pull/6))
- Improve input case support, VAPOIL/DISGAS systems, other fixes ([#8](https://github.com/sintefmath/JutulDarcy.jl/pull/8))

## [v0.1.5] - 2022-11-04

### Updates

- **Improved CSR support**
- **Performance improvements**

## [v0.1.4] - 2022-10-16

No specific changes documented.

## [v0.1.3] - 2022-09-29

### New Features

- **Region support for saturation tables**
- **Rewritten and documented main exported interfaces**
- **Support for choosing update frequency of AMG pressure subsystem**

## [v0.1.2] - 2022-09-25

### Improvements

- **Large improvements to compile time**
- **Added precompilation**
- **More tests**

## [v0.1.1] - 2022-09-20

### Merged Pull Requests

- Csr ([#1](https://github.com/sintefmath/JutulDarcy.jl/pull/1))
- Generalized AD ([#2](https://github.com/sintefmath/JutulDarcy.jl/pull/2))
- Refactor parameter logic + other fixes ([#3](https://github.com/sintefmath/JutulDarcy.jl/pull/3))
- Rewrite flux to support sensitivities ([#4](https://github.com/sintefmath/JutulDarcy.jl/pull/4))

### Closed Issues

- TagBot trigger issue ([#5](https://github.com/sintefmath/JutulDarcy.jl/pull/5))

