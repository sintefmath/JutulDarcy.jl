# Python interface

PyJutulDarcy provides a Python interface to JutulDarcy and manages the Julia installation. See the [PyJutulDarcy documentation](https://sintefmath.github.io/PyJutulDarcy/) for requirements and setup.

## Installation

Installation is done using `pip`. This will take care of all dependencies, and install a local Julia installation if it is not already present.

```shell
pip install jutuldarcy
```

## Usage

You can run simulations and convert the results to NumPy arrays. The [PyJutulDarcy user guide](https://sintefmath.github.io/PyJutulDarcy/guide/) covers indexing, results, Julia interoperation, execution backends, and gradients. To set up a simulation in Python, follow the [worked model-building example](https://sintefmath.github.io/PyJutulDarcy/examples/standalone/). See [Input formats](@ref) and [Supported physical systems](@ref) for details on the simulator.

```python
import jutuldarcy as jd
# Load SPE9 dataset to disk
pth = jd.test_file_path("SPE9", "SPE9.DATA")
# Simulate the model and convert to Python dicts
res = jd.simulate_data_file(pth, convert = True)
# Get field quantities and plot
import matplotlib.pyplot as plt
fopr = res["FIELD"]["FOPR"]
days = res["DAYS"]
plt.plot(days, fopr)
plt.ylabel("Field oil production")
plt.xlabel("Days")
plt.show()
```
