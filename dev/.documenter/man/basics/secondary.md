
# Secondary variables (properties) {#Secondary-variables-properties}

## Fluid systems {#Fluid-systems}

### General {#General}

#### Relative permeabilities {#Relative-permeabilities}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.BrooksCoreyRelativePermeabilities' href='#JutulDarcy.BrooksCoreyRelativePermeabilities'><span class="jlbinding">JutulDarcy.BrooksCoreyRelativePermeabilities</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BrooksCoreyRelativePermeabilities(
    sys_or_nph::Union{MultiPhaseSystem, Integer},
    exponents = 1.0,
    residuals = 0.0,
    endpoints = 1.0
)
```


Secondary variable that implements the family of Brooks-Corey relative permeability functions. This is a simple analytical expression for relative permeabilities that has a limited number of parameters:

$K(S) = K_{max} * ((S - S_r)/(1 - S_r^{tot}))^N$

**Fields**
- `exponents`: Exponents for each phase
  
- `residuals`: Residual saturations for each phase
  
- `endpoints`: Maximum relative permeability for each phase
  
- `residual_total`: Total residual saturation over all phases
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/relperm/simple.jl#L91-L107" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.RelativePermeabilities' href='#JutulDarcy.RelativePermeabilities'><span class="jlbinding">JutulDarcy.RelativePermeabilities</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
RelativePermeabilities((kr1, kr2, ...))
```


A simple relative permeability implementation. Assumes that each phase has a relative permeability on the form:

$K_{r,phase} = F(S_{phase})$

Supports multiple fluid regions through the `regions` keyword.

**Examples**

Single region:

```
kr1 = S -> S^2
kr2 = S -> S^3

kr = RelativePermeabilities((kr1, kr2))
```


Two regions:

```
kr1_reg1 = S -> S^2
kr2_reg1 = S -> S^3

kr1_reg2 = S -> S^3
kr2_reg2 = S -> S^4

regions # should be a vector with one entry that is 1 or 2 for each cell in the domain

kr = RelativePermeabilities(((kr1_reg1, kr2_reg1), (kr1_reg2, kr2_reg2)), regions = regions)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/relperm/simple.jl#L3-L33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.TabulatedSimpleRelativePermeabilities' href='#JutulDarcy.TabulatedSimpleRelativePermeabilities'><span class="jlbinding">JutulDarcy.TabulatedSimpleRelativePermeabilities</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TabulatedSimpleRelativePermeabilities(s::AbstractVector, kr::AbstractVector; regions::Union{AbstractVector, Nothing} = nothing, kwarg...)
```


Interpolated multiphase relative permeabilities that assumes that the relative permeability of each phase depends only on the phase saturation of that phase.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/relperm/simple.jl#L132-L137" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.ReservoirRelativePermeabilities' href='#JutulDarcy.ReservoirRelativePermeabilities'><span class="jlbinding">JutulDarcy.ReservoirRelativePermeabilities</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ReservoirRelativePermeabilities(
    w = nothing,
    g = nothing,
    ow = nothing,
    og = nothing,
    scaling = NoKrScale(),
    regions = nothing
    hysteresis_w = NoHysteresis(),
    hysteresis_ow = NoHysteresis(),
    hysteresis_og = NoHysteresis(),
    hysteresis_g = NoHysteresis(),
    hysteresis_s_threshold = 0.0,
    hysteresis_s_eps = 1e-10
)
```


Relative permeability with advanced features for reservoir simulation. Includes features like rel. perm. endpoint scaling, connate water adjustment and separate phase pair relative permeabilites for the oil phase.

**Fields**
- `krw`
  
- `krow`
  
- `krog`
  
- `krg`
  
- `regions`
  
- `phases`
  
- `hysteresis_w`
  
- `hysteresis_ow`
  
- `hysteresis_og`
  
- `hysteresis_g`
  
- `scaling`
  
- `hysteresis_s_threshold`
  
- `hysteresis_s_eps`
  

**Examples**

```julia
s = collect(range(0, 1, 100))
krw = PhaseRelativePermeability(s, s)
krog = PhaseRelativePermeability(s, s.^3)
kr_def = ReservoirRelativePermeabilities(krw = krw, krog = krog)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/relperm/advanced.jl#L31-L62" target="_blank" rel="noreferrer">source</a></Badge>

</details>


The `ReservoirRelativePermeabilities` type also supports hysteresis for either phase.
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.NoHysteresis' href='#JutulDarcy.NoHysteresis'><span class="jlbinding">JutulDarcy.NoHysteresis</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
NoHysteresis()
```


Type to indicate that no hysteresis is active, and the drainage curve will always be used.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/relperm/hysteresis.jl#L34-L39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.CarlsonHysteresis' href='#JutulDarcy.CarlsonHysteresis'><span class="jlbinding">JutulDarcy.CarlsonHysteresis</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CarlsonHysteresis()
```


Carlson&#39;s hysteresis model.

Note that this model requires an intersection between drainage and imbibition relative permeability that comes at some cost during simulation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/relperm/hysteresis.jl#L17-L24" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.KilloughHysteresis' href='#JutulDarcy.KilloughHysteresis'><span class="jlbinding">JutulDarcy.KilloughHysteresis</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
KilloughHysteresis(tol = 0.1, s_min = 0.0)
```


Killough hysteresis model. `tol` is a parameter for numerical stability and `s_min` a minimum threshold for when hysteresis is activated. Larger values for both of these parameters reduce numerical difficulties.

Reference: https://doi.org/10.2118/5106-PA


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/relperm/hysteresis.jl#L3-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.JargonHysteresis' href='#JutulDarcy.JargonHysteresis'><span class="jlbinding">JutulDarcy.JargonHysteresis</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
JargonHysteresis()
```


Jargon&#39;s hystersis model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/relperm/hysteresis.jl#L27-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.ImbibitionOnlyHysteresis' href='#JutulDarcy.ImbibitionOnlyHysteresis'><span class="jlbinding">JutulDarcy.ImbibitionOnlyHysteresis</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ImbibitionOnlyHysteresis()
```


Type to indicate that the hysteresis does not make use of the drainage curve, and the imbibition curve will always be used.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/relperm/hysteresis.jl#L42-L47" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.PhaseRelativePermeability' href='#JutulDarcy.PhaseRelativePermeability'><span class="jlbinding">JutulDarcy.PhaseRelativePermeability</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PhaseRelativePermeability(s, k; label = :w, connate = s[1], epsilon = 1e-16)
```


Type that stores a sorted phase relative permeability table (given as vectors of equal length `s` and `k`):

$K_r = K(S)$

Optionally, a label for the phase, the connate saturation and a small epsilon value used to avoid extrapolation can be specified. The return type holds both the table, the phase context, the autodetected critical and maximum relative permeability values and can be passed saturation values to evaluate the underlying function:

```julia
s = range(0, 1, 50)
k = s.^2
kr = PhaseRelativePermeability(s, k)
round(kr(0.5), digits = 2)

# output

0.25
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/types.jl#L275-L299" target="_blank" rel="noreferrer">source</a></Badge>

</details>


#### Phase viscosities {#Phase-viscosities}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.DeckPhaseViscosities' href='#JutulDarcy.DeckPhaseViscosities'><span class="jlbinding">JutulDarcy.DeckPhaseViscosities</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DeckPhaseViscosities(pvt, regions = nothing)
```


Secondary variable used to evaluate viscosities when a case is generated from a input file. Typically not instantiated in user scripts.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/deck_types.jl#L5-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.ConstMuBTable' href='#JutulDarcy.ConstMuBTable'><span class="jlbinding">JutulDarcy.ConstMuBTable</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConstMuBTable(pvtw::M) where M<:AbstractVector
```


Create a constant viscosity and formation-volume-factor table from a vector. Typical usage is to wrap a PVTW type table generated from external software.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/deck_types.jl#L159-L164" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.MuBTable' href='#JutulDarcy.MuBTable'><span class="jlbinding">JutulDarcy.MuBTable</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MuBTable(pvt, regions = nothing)
```


Table used to evaluate viscosities and shrinkage factors when a case is generated from a input file. Typically used to wrap tables (e.g. PVDG, PVDO) for use in simulation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/deck_types.jl#L77-L83" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.PhaseMassDensities' href='#JutulDarcy.PhaseMassDensities'><span class="jlbinding">JutulDarcy.PhaseMassDensities</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type representing the evaluation of mass density of each phase (i.e. units of mass per units of volume, for each cell in the model domain.)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/pvt.jl#L1-L4" target="_blank" rel="noreferrer">source</a></Badge>

</details>


#### Phase densities {#Phase-densities}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.DeckPhaseMassDensities' href='#JutulDarcy.DeckPhaseMassDensities'><span class="jlbinding">JutulDarcy.DeckPhaseMassDensities</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DeckPhaseMassDensities(pvt, regions = nothing)
```


Secondary variable used to evaluate densities when a case is generated from a input file. Typically not instantiated in user scripts.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/deck_types.jl#L29-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>


#### Shrinkage factors {#Shrinkage-factors}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.DeckShrinkageFactors' href='#JutulDarcy.DeckShrinkageFactors'><span class="jlbinding">JutulDarcy.DeckShrinkageFactors</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



DeckShrinkageFactors(pvt, regions = nothing)

Secondary variable used to evaluate shrinkage factors when a case is generated from a input file. Typically not instantiated in user scripts.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/deck_types.jl#L53-L58" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.ConstantCompressibilityDensities' href='#JutulDarcy.ConstantCompressibilityDensities'><span class="jlbinding">JutulDarcy.ConstantCompressibilityDensities</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConstantCompressibilityDensities(
    sys_or_nph::Union{MultiPhaseSystem, Integer},
    reference_pressure = 1.0,
    reference_density = 0.0,
    compressibility = 1.0
)
```


Secondary variable that implements a constant compressibility relationship for density. Given the reference pressure, compressibility and density at the reference pressure, each phase density can be computed as:

$ρ(S) = ρ_{ref} e^{(p - p_{ref})c}$

The constructor can take in either one value per phase or a single value for all phases for the reference pressure, compressibility and density at reference conditions.

**Fields**
- `reference_pressure`: Reference pressure for each phase (where the reference densities are given)
  
- `reference_densities`: Densities at the reference point
  
- `compressibility`: Compressibility factor used when expanding around reference pressure, typically between 1e-3 and 1e-10
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/pvt.jl#L8-L28" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Black-oil flow {#Black-oil-flow}

### Compositional flow {#Compositional-flow}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.PhaseMassFractions' href='#JutulDarcy.PhaseMassFractions'><span class="jlbinding">JutulDarcy.PhaseMassFractions</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PhaseMassFractions(:liquid)
```


Variable that defines the component mass fractions in a specific phase.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/multicomponent/variables/others.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.KValueWrapper' href='#JutulDarcy.KValueWrapper'><span class="jlbinding">JutulDarcy.KValueWrapper</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
KValueWrapper(K; dependence::Symbol = :pT)
```


Create a wrapper for a K-value interpolator to be used with K-value flash.

The main purpose of this wrapper is to transform the general flash cond NamedTuple into the right arguments for multi-linear interpolation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/multicomponent/utils.jl#L30-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Wells {#Wells}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.TotalMass' href='#JutulDarcy.TotalMass'><span class="jlbinding">JutulDarcy.TotalMass</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TotalMasses()
```


Variable that defines total mass of all components in each cell of the domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/multiphase.jl#L173-L177" target="_blank" rel="noreferrer">source</a></Badge>

</details>

