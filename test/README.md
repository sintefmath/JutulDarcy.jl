# Running tests

Tests use [Julia Test Items](https://julia-testitems.org/). CI uses
`julia-actions/julia-run-testitems@v2`; `Pkg.test()` remains supported.

With `juliati` installed, run a topic or combine topics:

```sh
juliati --filter ':wells in tags'
juliati --filter ':wells in tags && :adjoint in tags' --max-workers 1
juliati --filter ':thermal in tags && :backends in tags' --max-workers 1
juliati --filter ':ci_flow in tags' --max-workers 1 --check-bounds yes
```

Using `Pkg.test()` instead:

```julia
using Pkg
Pkg.test() # everything
ENV["JUTULDARCY_TEST_TAGS"] = "wells,adjoint" # require BOTH thematic tags
Pkg.test()
delete!(ENV, "JUTULDARCY_TEST_TAGS")
ENV["JUTULDARCY_TEST_GROUP"] = "flow"
Pkg.test()
delete!(ENV, "JUTULDARCY_TEST_GROUP")
```

Thematic tags are independent of the `:ci_*` scheduling tags. They describe the
coverage of an entire test item, so selecting `:wells` can also run other physics
covered by that item.

| Area | Tags |
| --- | --- |
| Physics | `:flow`, `:singlephase`, `:immiscible`, `:blackoil`, `:compositional`, `:thermal`, `:co2`, `:tracers` |
| Wells | `:wells`, `:controls`, `:groups`, `:multisegment` |
| Properties | `:properties`, `:systems`, `:relperm`, `:capillary`, `:tabulation`, `:flash`, `:diffusion` |
| Sensitivities | `:adjoint`, `:parameters`, `:derivatives`, `:optimization`, `:history_matching`, `:initialization`, `:equilibration` |
| Numerics | `:solvers`, `:cpr`, `:nldd`, `:sequential`, `:discretization`, `:flux`, `:sparsity`, `:sources` |
| Backends | `:backends`, `:kernel_abstractions`, `:gpu`, `:float32`, `:scalarization`, `:parallel`, `:mpi` |
| Geometry and input | `:mesh`, `:input`, `:fractures`, `:dfm`, `:utilities` |

## CI groups

Each item has exactly one `:ci_*` tag. The `CI group coverage` item checks that new
items are assigned to a group; `Pkg.test()` also validates tags while filtering.
Keep the matrix in `.github/workflows/CI.yml` and `ci_groups.jl` in sync.

| Group | Shared compilation / subject | Supplied serial timings |
| --- | --- | ---: |
| `core` | Mesh, utilities, properties, well merging | 12m29s |
| `flow` | Reservoir simulations, input decks, discretizations | 9m18s |
| `wells` | Coupled wells and facility controls | 12m38s |
| `backends` | KernelAbstractions and Float32 simulation | 12m57s |
| `solvers` | CPR, NLDD and fractured flow | 12m46s |
| `partitioned` | PartitionedArrays and MPI solver extensions | 8m01s |
| `thermal` | Thermal simulation and tabulated properties | 9m59s |
| `compositional` | Compositional/CO2 physics and thermal/tracer backends | 12m41s |
| `adjoint_parameters` | Parameter optimization, state and equilibrium gradients | 11m57s |
| `adjoint_controls` | Inner adjoints and force/control gradients | 11m35s |

These are sums from the supplied 114-minute serial run, not measured times for
fresh CI workers. The target is 15 minutes per group including startup. Each job
uses one test process so related items share compiled code. Results JSON artifacts
and the aggregate report provide timings for tuning after the first CI run.
One preparation job resolves the package manifest and warms the full test sandbox
by running only the `:infrastructure` item with the same coverage and forced bounds
checks as the groups. It uses a portable CPU target, then all groups restore its
exact depot cache and package manifest. Installation and test workers consistently
use the `release` channel; the resolved version is recorded in the cache key. The
runner's own toolkit depot is separate from this test depot. Group jobs do not
create ten redundant copies of that cache.

A cold preparation can exceed the 15-minute target; its timeout is 45 minutes.
The test groups retain a 30-minute timeout while their actual runtimes are measured.
Reactivating the prepared environment has a separate 10-minute limit so a startup
failure can produce test results and diagnostics before the job timeout kills it.

The multimodel sensitivity item was split at its existing `Inner adjoint tests`
testset so its original 16-minute runtime does not force one group over target.
Existing failures and broken assertions remain in place.

The workflow peaks at 12 concurrent jobs: 10 test groups, MPI and documentation.
The results report runs after the test matrix completes.
