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
```

### Imposing limits on wells (multiple constraints)

## Well forces

### Changing perforations

```@docs
PerforationMask
```

### Other forces

Can use [`SourceTerm`](@ref) or [`FlowBoundaryCondition`](@ref)
