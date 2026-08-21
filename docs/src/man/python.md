# Python interface

JutulDarcy has a Python interface that automatically manages the Julia installation. For more details, see the [PyJutulDarcy GitHub page](https://github.com/sintefmath/PyJutulDarcy).

## Installation

Installation is done using `pip`. This will take care of all dependencies, and install a local Julia installation if it is not already present.

```shell
pip install jutuldarcy
```

## Usage

The package contains a high level interface that allows you to run simulations and convert the results to `numpy` arrays. In addition, the package uses [PythonCall](https://github.com/JuliaPy/PythonCall.jl) under the hood where it is possible to access functions that do not yet have a public API.

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
