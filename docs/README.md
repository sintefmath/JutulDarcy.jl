# Building documentation

Instantiate the docs environment, then build and execute the Literate pages before
assembling the site:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=docs docs/make.jl --examples-only
julia --project=docs docs/make.jl --docs-only
```

On headless Linux, prefix the build commands with `xvfb-run -a`. Running
`julia --project=docs docs/make.jl` performs both stages automatically.

`JUTULDARCY_DOCS_EXAMPLE_WORKERS` controls concurrency (default: 2). Julia tasks
launch isolated Julia processes with one Julia thread and one OpenBLAS thread
each. Each process builds complete example categories to share compilation.
Isolation is necessary because working directories, output capture and GLMakie
state are shared within a process. Two workers bound the additional memory use;
increase the setting on machines with enough memory.

Literate executes the scripts and saves static Markdown and figure assets in
`docs/src/examples/`. The final Documenter pass consumes those results without
executing the example scripts again. Inline examples in the manual still run in
Documenter, and notebooks are generated as before. Worker logs are in
`docs/example_logs/` and are uploaded by CI even when an example fails. Worker
output is also streamed into the CI log, with a start message for every example.
Failed workers fail the build. Successful pages have source fingerprints;
`--docs-only` rejects missing or stale pages.
Embedded HTML PNG/JPEG images are saved as separate assets. Text output is bounded
to 80 lines and 160 columns because Literate does not apply Documenter's limited
display context to large returned tables; truncated output is labeled explicitly.

The existing example selection and draft environment settings remain supported.
For interactive builds, set `JUTULDARCY_DOCS_AUTOBUILD=0` before including
`docs/make.jl`, then call `build_jutul_darcy_docs` with the desired keywords.
Run the prebuild regression checks with:

```sh
julia --project=docs docs/test/runtests.jl
```

## Further CI splitting

The two stages also allow an artifact-based matrix later: build disjoint sets of
categories on separate runners, upload their generated Markdown, figures and
fingerprints, then download and merge them before `--docs-only`. To build a subset
without assembling a partial site:

```julia
ENV["JUTULDARCY_DOCS_AUTOBUILD"] = "0"
include("docs/make.jl")
categories = ("introduction", "properties")
examples = [(category, name) for (category, names) in get_example_paths()
            for name in names if category in categories]
prebuild_examples(examples; workers = 2)
```

Use artifacts from the same workflow run and share the resolved docs Manifest
between producers and the final job. Upload `docs/src/examples/` after successful
execution; generate the overview in the final job. Avoid reusing executed pages
based only on the example script: simulator and dependency changes affect results.
Fingerprints include the simulator source, build scripts, Julia version and docs
Manifest when present.

A four-job example matrix would bring the workflow peak to 15 jobs (10 tests,
MPI and four example producers); the final docs job would follow the producers.
Measure the two-worker build first to decide whether the extra environment setup
and artifact transfers are worthwhile under the 20-job concurrency limit.

After successful deployment, CI creates or updates one preview-link comment on
same-repository PRs. Fork PRs retain the existing read-only-token deployment
restrictions; the workflow does not grant their code additional permissions.
