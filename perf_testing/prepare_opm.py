"""Export the same two layer-refined benchmarks through OPM's deck parser.

Run with the system Python that provides opm.io and numpy. This deliberately
supports only these benchmark decks, following refine_layers.jl. Source data
are read only; all generated decks and manifests belong under results/.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from opm.io.parser import Parser, ParseContext, action


def values(record):
    result = []
    for item in record:
        for index in range(len(item)):
            if getattr(item, "__defaulted")(index):
                value = None
            elif item.is_int():
                value = item.get_int(index)
            elif item.is_double():
                value = item.get_raw(index)
            elif item.is_string():
                value = item.get_str(index)
            else:
                raise ValueError(f"Unsupported transformed item: {item.name()}")
            result.append(value)
    return result


def token(value):
    if value is None:
        return "1*"
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    return str(value) if isinstance(value, int) else format(value, ".17g")


def records_text(name, rows, terminated=True):
    return name + "\n" + "".join(" ".join(map(token, row)) + " /\n" for row in rows) + ("/\n" if terminated else "")


def array_text(name, array):
    flat = np.asarray(array).ravel(order="F")
    return name + "\n" + "\n".join(" ".join(map(token, flat[i:i+8].tolist())) for i in range(0, len(flat), 8)) + "\n/\n"


def table_text(keyword, width):
    text = keyword.name + "\n"
    for record in keyword:
        if keyword.name in ("PVTO", "PVTG"):
            if len(record[1]) == 0:
                text += "/\n"
                continue
            row = values(record)
            text += " ".join(map(token, row[:4])) + "\n"
            row = row[4:]
        else:
            row = values(record)
        assert len(row) % width == 0
        text += "".join(" ".join(map(token, row[i:i+width])) + "\n" for i in range(0, len(row), width)) + "/\n"
    # Empty records separate PVT regions; the final region terminator is not
    # exposed as a record by OPM's Python iterator.
    return text + ("/\n" if keyword.name in ("PVTO", "PVTG") else "")


def prepare(source, target, factor, fraction):
    context = ParseContext()
    # Legacy Olympus has COORDSYS without NUMRES and a redundant EDIT slash.
    # Jutul ignores COORDSYS. No numerical property is omitted here.
    context.ignore_keyword("COORDSYS")
    context.update("PARSE_MISSING_DIMS_KEYWORD", action.warn)
    context.update("PARSE_RANDOM_SLASH", action.warn)
    parser = Parser()
    deck = parser.parse(str(source), context)
    dims = tuple(values(deck["DIMENS"][0])[:3])
    nx, ny, nz = dims
    cells = nx * ny * nz
    newdims = (nx, ny, nz * factor)
    forbidden = {"NNC", "TRANX", "TRANY", "TRANZ", "AQUNUM", "AQUCON", "COMPSEGS", "COMPDATMD", "WPIMULT", "TSTEP"}
    assert not any(kw.name in forbidden for kw in deck)
    reports = sum(len(kw) for kw in deck if kw.name == "DATES")
    selected_reports = math.ceil(fraction * reports)
    report = 0
    arrays = {}
    table_widths = dict(SWOF=4, SGOF=4, SLGOF=4, SWFN=3, SGFN=3,
        SOF3=3, SOF2=2, PVDO=3, PVDG=3, PVTO=3, PVTG=3,
        RSVD=2, RVVD=2, PBVD=2, PDVD=2, TEMPVD=2)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w") as output:
        output.write("-- Layer refinement matching perf_testing/refine_layers.jl\n")
        for kw in deck:
            name = kw.name
            if name == "GRID" and "EQLDIMS" not in deck and "EQUIL" in deck:
                output.write(records_text("EQLDIMS", [[len(deck["EQUIL"])]], False))
            if name in ("END", "COORDSYS"):
                continue
            if name == "DATES":
                assert len(kw) == 1
                output.write(str(kw))
                report += 1
                if report == selected_reports:
                    break
                continue
            if name == "GRID":
                output.write(str(kw))
                if "MINPV" not in deck and "MINPORV" not in deck:
                    # GeoEnergyIO inserts this metric-deck default before the
                    # Julia refinement helper divides it among child layers.
                    assert "METRIC" in deck
                    output.write(records_text("MINPV", [[1e-6/factor]], False))
            elif name == "WELLDIMS":
                row = values(kw[0])
                if row[1] is not None and row[1] > 0:
                    row[1] *= factor
                output.write(records_text(name, [row], False))
            elif name in ("DIMENS", "SPECGRID"):
                row = values(kw[0])
                row[:3] = newdims
                output.write(records_text(name, [row], False))
            elif name in table_widths:
                output.write(table_text(kw, table_widths[name]))
            elif name == "ZCORN":
                old = kw.get_raw_array().reshape((2*nx, 2*ny, 2*nz), order="F")
                refined = np.empty((2*nx, 2*ny, 2*nz*factor))
                for k in range(nz):
                    top, bottom = old[:, :, 2*k], old[:, :, 2*k+1]
                    for sub in range(factor):
                        child = k*factor + sub
                        refined[:, :, 2*child] = top + (sub/factor)*(bottom-top)
                        refined[:, :, 2*child+1] = top + ((sub+1)/factor)*(bottom-top)
                output.write(array_text(name, refined))
                arrays[name] = refined
            elif len(kw) == 1 and len(kw[0]) == 1 and len(kw[0][0]) == cells:
                item = kw[0][0]
                old = kw.get_int_array() if item.is_int() else kw.get_raw_array()
                old = old.reshape(dims, order="F")
                refined = np.repeat(old, factor, axis=2)
                if name == "PORV":
                    refined = refined / factor
                elif name in ("MULTZ", "MULTZ-"):
                    refined.fill(1)
                    offset = factor-1 if name == "MULTZ" else 0
                    refined[:, :, offset::factor] = old
                output.write(array_text(name, refined))
                arrays[name] = refined
            elif name == "MINPV":
                output.write(records_text(name, [[kw[0][0].get_raw(0)/factor]], False))
            elif name in ("FAULTS", "COMPDAT", "WELOPEN", "BOX", "EQUALS", "MULTIPLY", "COPY", "ADD"):
                rows = []
                for record in kw:
                    row = values(record)
                    if name == "FAULTS":
                        low, high = row[5:7]
                        if row[7].startswith(("Z", "K")):
                            for k in range(low, high+1):
                                child = (k-1)*factor+1 if row[7].endswith("-") else k*factor
                                rows.append(row[:5]+[child, child]+row[7:])
                            continue
                        row[5:7] = [(low-1)*factor+1, high*factor]
                    elif name == "COMPDAT":
                        row[3:5] = [(row[3]-1)*factor+1, row[4]*factor]
                        for index in (7, 9):
                            if row[index] is not None and row[index] > 0:
                                row[index] /= factor
                    elif name == "WELOPEN":
                        assert all(value is None or value <= 0 for value in row[5:7])
                        if row[4] is not None and row[4] > 0:
                            for sub in range(factor):
                                rows.append(row[:4]+[(row[4]-1)*factor+sub+1]+row[5:])
                            continue
                    else:
                        start = 4 if name == "BOX" else 6
                        if name != "BOX" and row[0] in ("MULTZ", "MULTZ-"):
                            assert name == "EQUALS", "Unvalidated interface multiplier operation"
                            low = 1 if row[start] is None else row[start]
                            high = nz if row[start+1] is None else row[start+1]
                            for k in range(low, high+1):
                                child = k*factor if row[0] == "MULTZ" else (k-1)*factor+1
                                rows.append(row[:start]+[child, child]+row[start+2:])
                            continue
                        if row[start] is not None:
                            row[start] = (row[start]-1)*factor+1
                        if row[start+1] is not None:
                            row[start+1] *= factor
                    rows.append(row)
                output.write(records_text(name, rows, name != "BOX"))
            else:
                output.write(str(kw))
        output.write("END\n")
    assert report == selected_reports
    # Round-trip arrays through the same parser before permitting a Flow run.
    check = parser.parse(str(target))
    roundtrip_errors = {}
    for name, expected in arrays.items():
        actual = check[name].get_int_array() if np.issubdtype(expected.dtype, np.integer) else check[name].get_raw_array()
        expected = expected.ravel(order="F")
        np.testing.assert_allclose(actual, expected, rtol=5e-16, atol=0)
        nonzero = expected != 0
        roundtrip_errors[name] = float(np.max(np.abs((actual[nonzero]-expected[nonzero])/expected[nonzero]), initial=0))
    manifest = dict(source=str(source), target=str(target), factor=factor,
        dimensions=list(map(int, newdims)), reports=selected_reports,
        source_reports=reports, fraction=fraction,
        active_cells=int(np.count_nonzero(arrays["ACTNUM"])),
        array_roundtrip_relative_errors=roundtrip_errors,
        parser_compatibility=["COORDSYS omitted", "redundant slash accepted", "missing EQLDIMS made explicit"])
    target.with_suffix(".json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest), flush=True)


if __name__ == "__main__":
    arguments = argparse.ArgumentParser()
    arguments.add_argument("source", type=Path)
    arguments.add_argument("target", type=Path)
    arguments.add_argument("factor", type=int)
    arguments.add_argument("--fraction", type=float, default=0.2)
    args = arguments.parse_args()
    assert args.factor >= 1 and 0 < args.fraction <= 1
    prepare(args.source, args.target, args.factor, args.fraction)
