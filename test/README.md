# Running tests

`Pkg.test()` runs the full suite by default. Set `JUTULDARCY_TEST_GROUP` to a
group from the CI matrix to run that shard. Set `JUTULDARCY_TEST_TAGS` to a
comma-separated list to run items that have every listed tag. For example:

```sh
JUTULDARCY_TEST_GROUP=thermal_input julia --project -e 'using Pkg; Pkg.test()'
JUTULDARCY_TEST_TAGS=wells,adjoint julia --project -e 'using Pkg; Pkg.test()'
```

Editors and other Julia Test Items runners can filter by the same tags, such
as `:tracers` or `:ci_thermal_input`.

Give each new `@testitem` exactly one `:ci_*` tag from `ci_groups.jl` and any
thematic tags that help people find it. The CI matrix uses these groups to
run related tests in the same Julia process.
