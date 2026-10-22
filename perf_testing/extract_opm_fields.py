"""Extract OPM's final reservoir fields with explicit global-cell identities."""
import json
import sys
from pathlib import Path
import numpy as np
from opm.io.ecl import ERst, EclFile

restart_path, prefix = map(Path, sys.argv[1:3])
rst = ERst(str(restart_path))
report = max(rst.report_steps)
assert int(rst['INTEHEAD', report][2]) == 1, 'This comparison expects metric-deck output'
pressure_source = rst['PRESSURE', report]
water_source = rst['SWAT', report]
gas_source = rst['SGAS', report] if ('SGAS', report) in rst else np.zeros_like(water_source)
pressure = np.asarray(pressure_source, dtype=np.float64)*1e5
water, gas = map(lambda a: np.asarray(a, dtype=np.float64), (water_source, gas_source))
saturations = np.vstack((water, 1-water-gas, gas))
grid = EclFile(str(restart_path.with_suffix('.EGRID')))
cells = np.flatnonzero(grid['ACTNUM']).astype('<i8')+1
assert len(cells) == len(pressure), (len(cells), len(pressure))
assert np.all(np.isfinite(pressure)) and np.all(pressure > 0)
assert np.all(np.isfinite(saturations))
prefix.parent.mkdir(parents=True, exist_ok=True)
pressure.astype('<f8').tofile(str(prefix)+'.pressure.bin')
saturations.astype('<f8').ravel(order='F').tofile(str(prefix)+'.saturations.bin')
cells.tofile(str(prefix)+'.cells.bin')
metadata = dict(cells=len(cells), report=int(report), elapsed_days=float(rst['DOUBHEAD',report][0]), pressure_unit='Pa',
    phase_names=['Aqueous','Liquid','Vapor'],
    pressure_source_type=str(pressure_source.dtype),
    water_source_type=str(water_source.dtype), gas_source_type=str(gas_source.dtype))
Path(str(prefix)+'.fields.toml').write_text(''.join(f'{key} = {json.dumps(value)}\n' for key,value in metadata.items()))
print(json.dumps(metadata), flush=True)
