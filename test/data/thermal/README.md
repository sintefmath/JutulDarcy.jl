These self-contained DATA cases are reduced analogues of the thermal/tracer
examples used to exercise the Eclipse input workflow:

| Deck | Grid | Coverage |
| --- | --- | --- |
| `HEAT_SMALL.DATA` | 5 x 5 x 2 | Three-phase dissolved-gas black oil, cold-water injection, SPECROCK followed by VISCREF, all three temperature-dependent viscosity tables, nonuniform geometry |
| `THERMAL_TRACER_SMALL.DATA` | 5 x 7 x 3 | Oil/water five-spot, RTEMPVD initialization, fluid and rock heat capacities, passive water tracer, temperature/tracer updates without a repeated rate control |
| `TRACER_SMALL.DATA` | 5 x 5 x 1 | Water-only quarter-five-spot, direct pressure initialization without SWAT, omitted vertical permeability, tracer injection and transport |

These cases require the companion GeoEnergyIO.jl parser changes on
[`thermal-keyword-support`](https://github.com/moyner/GeoEnergyIO.jl/tree/thermal-keyword-support).
Develop that checkout in the test environment when testing this branch.

Run them with `setup_case_from_data_file` and `simulate_reservoir`. The high
pressure black-oil case requests a 1500 bar BHP change at startup. Use
`setup_case_from_data_file(path; dp_max_abs_well = Inf)` so the default 50 bar
Newton update cap does not prevent the well from reaching its target within
the iteration limit. Timesteps still adapt to nonlinear convergence.
The automated checks are in `test/thermal_data.jl` and are registered in
`test/runtests.jl`. They use tighter mass and tracer convergence tolerances,
and check tracer bounds weighted by carrier-phase saturation in nearly dry cells.
The original grids have 1200, 56100 and 225 cells, respectively. Rates and
schedules have been reduced for short regression runs; these are not reference
solutions for the full models.

The `DIFF` option in `TRACER_SMALL.DATA` deliberately checks the diagnostic for
unsupported numerical dispersion control. Transport uses JutulDarcy's existing
implicit upwind tracer scheme. It does not reproduce Eclipse's dispersion-control
algorithm. Tracers currently support the free carrier phase and default units;
partitioning, recycling and cumulative injection factors are rejected.

`TEMP` uses internal-energy transport with pressure work omitted, following
the [OPM temperature formulation](https://opm.github.io/opm-flow-editor-support/keywords/TEMP.html).
`THERMAL` retains pressure work. These tests check the implemented formulation;
they do not establish equivalence with either commercial simulator formulation.
Temperature-dependent heat capacities are integrated into internal energy,
so differentiating the energy with respect to temperature recovers the input
heat capacity.

Thermal units follow the [OPM keyword reference](https://opm.github.io/opm-flow-editor-support/keywords/SPECROCK.html):
SPECROCK is heat capacity per rock volume, while SPECHEAT is per fluid mass.
[THCONR](https://opm.github.io/opm-flow-editor-support/keywords/THCONR.html)
is the bulk rock-plus-fluid conductivity and is not scaled again by rock fraction.
