# Wells and controls

## Well setup routines

Wells can be set up using the convenience functions [`setup_well`](@ref) and [`setup_vertical_well`](@ref). These routines act on the output from [`reservoir_domain`](@ref) and can set up both types of wells. We recommend that you use these functions instead of manually calling the well constructors.

```@docs
JutulDarcy.WellGroup
```

## Types of wells

### Simple wells

```@docs
JutulDarcy.SimpleWell
```

#### Equations

### Multisegment wells

```@docs
WellDomain
MultiSegmentWell
JutulDarcy.SegmentWellBoreFrictionHB
JutulDarcy.SegmentValve
JutulDarcy.PotentialDropBalanceWell
MixedWellSegmentFlow
```

## Well controls and limits

### Types of well controls

```@docs
InjectorControl
ProducerControl
DisabledControl
```

```@docs
JutulDarcy.replace_target
JutulDarcy.default_limits
```

### Types of well targets

```@docs
BottomHolePressureTarget
SinglePhaseRateTarget
SurfaceLiquidRateTarget
SurfaceOilRateTarget
SurfaceGasRateTarget
SurfaceWaterRateTarget
TotalRateTarget
HistoricalReservoirVoidageTarget
ReservoirVoidageTarget
DisabledTarget
JutulDarcy.TotalReservoirRateTarget
TotalMassRateTarget
```

### Group rate targets

Pass `well_groups = Dict(:Group => [:W1, :W2])` to `setup_reservoir_model`
and apply the same group target to each member. Rates are positive for injection
and negative for production. Each group must contain only injectors or only
producers. Wells at individual limits contribute their current phase rates;
the remaining target is shared by wells still under group control.
Group targets are retained in the well limits under `group_rate`, `group_lrat`,
`group_orat`, `group_grat` or `group_wrat`, so wells can switch back to group
control when their allocated group rate becomes limiting.
Group rates, cumulative volumes and ratios are included in
`result.summary["VALUES"]["GROUP"]`. Use `plot_summary(result;
plots = ["Injection:GWIR", "Production:GWPR"])` to plot named groups.

`group_use_potentials = Dict(:Group => false)` gives equal weights instead of
the default potential weighting. This is stored in `WellGroup.use_potentials`.
See the example [Geothermal example with group controls](@ref).

```@docs
GroupTotalRateTarget
GroupSurfaceLiquidRateTarget
GroupSurfaceOilRateTarget
GroupSurfaceGasRateTarget
GroupSurfaceWaterRateTarget
```

### Implementation of well controls

```@docs
JutulDarcy.well_target_value
JutulDarcy.well_control_equation
JutulDarcy.check_well_limit
JutulDarcy.SurfaceComponentRates
```

### Well outputs

```@docs
JutulDarcy.print_well_result_table
```

### Imposing limits on wells (multiple constraints)

## Well forces

### Perforations and WI adjustments

```@docs
PerforationMask
JutulDarcy.Perforations
compute_peaceman_index
```

### Other forces

Can use [`SourceTerm`](@ref) or [`FlowBoundaryCondition`](@ref)
