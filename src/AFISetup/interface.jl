
"""
    case = setup_case_from_afi(afi_path)
    case = setup_case_from_afi(afi_dict)
    case = setup_case_from_afi(afi, reservoir = existing_reservoir)

Set up a `JutulCase` from an AFI input file. The input can be provided either as
a file path to the AFI file or as a dictionary representation of the AFI input
(from `read_afi_file` in GeoEnergyIO).

The parser requires that the mesh and mesh properties are either provided inline
(i.e. as direct keywords with numbers) or in the [RESQML
format](https://energistics.org/resqml-data-standards). GSG files are not
currently supported.

# Keyword arguments
- `step_limit`: An optional limit on the number of simulation steps to set up.
- `step_override`: An optional vector of time points (in seconds) to override the
  simulation steps defined in the AFI file. If provided, the simulation will use
  these time points as the time steps. If the time points are normalized (i.e. in the
  range [0, 1]), set `step_override_is_normalized` to true.
- `step_override_is_normalized`: A boolean indicating whether the `step_override`
  time points are normalized (i.e. in the range [0, 1]). Defaults to false.
- `reservoir`: An optional existing reservoir domain to use instead of setting
  up a new one from the AFI file.
"""
function setup_case_from_afi

end

function setup_case_from_afi(pth::AbstractString; verbose = false, kwarg...)
    x = read_afi_file(pth, convert = true, verbose = verbose)
    return setup_case_from_afi(AFIInputFile(x); kwarg...)
end

function setup_case_from_afi(x::AbstractDict; kwarg...)
    return setup_case_from_afi(AFIInputFile(x); kwarg...)
end

function setup_case_from_afi(afi::AFIInputFile;
        step_limit = missing,
        step_override = missing,
        step_override_is_normalized = false,
        kwarg...
    )
    model, prm = setup_reservoir_model(afi; extra_out = true, kwarg...)
    state0 = setup_reservoir_state(afi, model)
    dt, forces = setup_afi_schedule(afi, model,
        step_limit = step_limit,
        step_override = step_override,
        step_override_is_normalized = step_override_is_normalized,
    )
    date = first(keys(afi.setup["IX"]["STEPS"]))
    return Jutul.JutulCase(model, dt, forces, state0 = state0, parameters = prm, start_date = date)
end
