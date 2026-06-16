# Unit support

JutulDarcy uses consistent units. This typically means that all values must be
input in strict SI when you are defining a model. The advantage is that the code
itself is free from conversion factors, but this also means that you must
convert your units upfront when setting up a model.

## Input files

Input files (e.g. .DATA or .AFI files) are with the default settings automatically converted to SI when parsed. In this section of the manual, we highlight a few of the tools that can be used to handle units in your own scripts, as well as a few typical stumbling blocks.

## Comparison of unit systems

This section contains a brief comparison of the units in SI versus metric and field units for some key quantities:

| Quantity             | SI        | Metric     | Field                | Note                        |
|----------------------|-----------|------------|----------------------|------------------------------
| Length               | meter     | meter      | feet                 | Affects definition of volumes and rates
| Pressure             | pascal    | psi        | bar                  | Pascal is a very small unit (1 bar = 1e5 Pa)
| Time                 | second    | day        | day                  | Time is used in viscosity and permeability, which makes these values much smaller than normal
| Mass                 | kg        | kg         | pound                | Impacts densities of all kinds
| Reservoir volume     | meter$^3$ | meter$^3$  | stb                  | Volumes in in-situ conditions (varying pressure and temperature)
| Surface volume       | meter$^3$ | meter$^3$  | MScf (1000 ft$^3$) | Surface or standard conditions at some specified pressure and temperature. Field gas volume unit differs at surface and reservoir conditions
| Absolute temperature | Kelvin    | Kelvin     | Rankine              | Code uses absolute temperatures internally
| Relative temperature | Celsius   | Celsius    | Fahrenheit           | Important to convert using `convert_to_si` since there is both a shift and a scaling involved
| Permeability         | meter$^2$ | millidarcy | millidarcy           | The SI unit is *very* small ($m^2 \approx 10^{-12}$ Darcy)

From these values, additional units are defined. For instance the unit conversion factor for areas and volumes can be computed as $L^2$ and $L^3$, respectively where L is the unit for length. Similarly, the unit for mass density can be written as the unit for mass divided by the unit for volume and transmissibility is the product of the units for viscosity and reservoir volume divided by the product of the factors for time and pressure. More definitions used in conversions between unit systems can be found in the implementation of the [`GeoEnergyIO.InputParser.DeckUnitSystem type in GeoEnergyIO.jl`](https://github.com/sintefmath/GeoEnergyIO.jl/blob/main/src/InputParser/utils.jl).

### Typical pitfalls

Simulations typically output pressures, temperatures and surface rates. Surface rates are, for non-compositional models, just a scaling of the component mass rate done by dividing by the corresponding surface density. As SI rates are given per second, the values will be small. Pressures are given in Pascals, which will give large values. Permeability values are nominally small in magnitude for m$^2$, so it is important to avoid inputting millidarcy values directly. Wrong units is the most common reason for simulations not producing sensible results, or being impossible to simulate.

## Dealing with units in scripts

As mentioned, reading of input files will automatically convert data to the correct units for simulation, but care must be taken when you are writing your own code. `Jutul.jl` contains unit conversion factors to make it easier to write code. The main functions to use are [`Jutul.si_unit`](@ref), [`Jutul.convert_to_si`](@ref) and [`Jutul.convert_from_si`](@ref).

### Cheat sheet for units

Getting a named unit:

```@example
using Jutul
p = convert_to_si(120.0, :bar)
```

This can also be done with a `String` (useful if you are working in Python):

```@example
using Jutul
p = convert_to_si(120.0, "bar")
```

We can also use a magic `si` string to write this in a more compact form:

```@example
using Jutul
p = si"bar"
```

### Composite units

You can also extract individual units and compute conversion factors yourself:

```@example
using Jutul
day, stb = si_units(:day, :stb)
# convert to m^3/s:
rate = 100.0stb/day
```

The magic `si` string can also do this for you if you are writing a script:

```@example
using Jutul
rate = si"100stb/day"
```

### Internals

As there are no conversion factors internally in the code, you can in principle
use any consistent unit system. Some default scaling of variables assume that
the magnitude pressures and velocities roughly match that of strict SI (e.g.
Pascals and cubic meters per second). These scaling factors are primarily used
when iterative linear solvers are used. We recommend sticking to SI units.
