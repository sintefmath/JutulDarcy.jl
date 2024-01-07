# Wells and controls

## Types of wells

### Simple wells

```@docs
SimpleWell
```

#### Equations

### Multisegment wells

```@docs
MultiSegmentWell
```

### Well setup routines

[`setup_well`](@ref)
[`setup_vertical_well`](@ref)

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