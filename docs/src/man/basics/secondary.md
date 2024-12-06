# Secondary variables (properties)

## Fluid systems

### General

#### Relative permeabilities

```@docs
BrooksCoreyRelativePermeabilities
RelativePermeabilities
JutulDarcy.TabulatedSimpleRelativePermeabilities
JutulDarcy.ReservoirRelativePermeabilities
```

The `ReservoirRelativePermeabilities` type also supports hysteresis for either phase.

```@docs
JutulDarcy.NoHysteresis
JutulDarcy.CarlsonHysteresis
JutulDarcy.KilloughHysteresis
JutulDarcy.JargonHysteresis
JutulDarcy.ImbibitionOnlyHysteresis
```

```@docs
PhaseRelativePermeability
```

#### Phase viscosities

```@docs
DeckPhaseViscosities
JutulDarcy.ConstMuBTable
JutulDarcy.MuBTable
```

```@docs
PhaseMassDensities
```

#### Phase densities

```@docs
JutulDarcy.DeckPhaseMassDensities
```

#### Shrinkage factors

```@docs
DeckShrinkageFactors
```

```@docs
ConstantCompressibilityDensities
```

### Black-oil flow

```@docs
```

### Compositional flow

```@docs
PhaseMassFractions
JutulDarcy.KValueWrapper
```

## Wells

```@docs
TotalMass
```
