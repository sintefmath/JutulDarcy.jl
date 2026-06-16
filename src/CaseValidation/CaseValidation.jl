module CaseValidation
    using Jutul, JutulDarcy, Printf
    const VALDICT = Dict{Symbol, Tuple{Float64, Float64}}
    @enum SEVERITY OK = 0 WARNING = 1 ERROR = 2
    function validate(case::JutulCase;
            extrema_range = VALDICT(),
            mean_range = VALDICT(),
            info_level = 1,
            throw = false,
            eager = false
        )
        extrema = VALDICT()
        mean = VALDICT()
        default_units = Dict{Symbol, String}()
        custom_messages = Dict{Symbol, String}()

        function add_bounds(k, label::Union{String, Missing}, ext::Tuple{Float64, Float64}, mval::Union{Tuple{Float64, Float64}, Missing} = missing; msg = missing)
            if !ismissing(label)
                default_units[k] = label
            end
            extrema[k] = ext
            if !ismissing(mval)
                mean[k] = mval
            end
            if !ismissing(msg)
                custom_messages[k] = msg
            end
        end
        p0 = JutulDarcy.DEFAULT_MINIMUM_PRESSURE

        add_bounds(:permeability, "meters^2",
            (0.0, convert_to_si(10.0, :darcy)),
            (0.0, convert_to_si(5.0, :darcy)),
            msg = "Large values for permeability often indicates that that the values may have been input as millidarcy or darcy instead of m^2. Use `val = convert_to_si(val, \"millidarcy\")` to convert from darcy to m^2."
        )

        add_bounds(:porosity, missing,
            (1e-6, 1.0),
            (1e-6, 0.5),
            msg = "Porosity should be in the unit range. Zero porosity cells should be removed by passing `min_porevolume = 1e-6` when calling `reservoir_domain` during setup. Zero values will lead to singular systems."
        )

        add_bounds(:net_to_gross, missing,
            (1e-6, 1.0),
            msg = "Values should be in the unit range. Values larger than one are possible numerically, but may be unphysical. Zero net_to_gross cells should be removed by passing `min_porevolume = 1e-6` when calling `reservoir_domain` during setup. Zero values will lead to singular systems."
        )

        # Variables/parameters
        add_bounds(:Pressure, "Pascal", (p0, 5000*si_unit(:bar)), (p0, 1000*si_unit(:bar)), msg = "Pressures are given in Pascals. Large values may indicate bar, MPa or psi units. Use `val = convert_to_si(val, \"bar\")` to convert from bar to SI.")
        add_bounds(:Temperature, "Celsius", (0.0, convert_to_si(1000.0, :Celsius)))
        add_bounds(:Saturations, missing, (0.0, 1.0))
        add_bounds(:OverallMoleFractions, missing, (0.0, 1.0))
        add_bounds(:StaticFluidVolume, "meter^3", (0.01, 1e20))
        add_bounds(:FluidVolume, "meter^3", (0.01, 1e20))
        add_bounds(:WellIndices, missing, (0.00, 1e-6))

        # Well controls
        add_bounds(:rate, missing, (0.00, 100.0), msg = "Rates above 100 m^3/s or 100 kg/s are unlikely to be physical for a single well. Make sure that rates are not given per day.")
        add_bounds(:bhp, missing, (p0, 5000*si_unit(:bar)), msg = "Bottom-hole-pressures should be given in Pascals. Large values may indicate bar, MPa or psi units. Use `val = convert_to_si(val, \"bar\")` to convert from bar to SI.")

        for (k, v) in pairs(extrema_range)
            extrema[k] = v
        end
        for (k, v) in pairs(mean_range)
            mean[k] = v
        end
        result = (
            extrema = extrema,
            mean = mean,
            custom_messages = custom_messages,
            default_units = default_units,
            info_level = info_level,
            throw = throw,
            eager = eager,
            messages = Tuple{SEVERITY, String, String}[]
        )
        model = case.model
        print_progress(msg) = info_level > 0 && print_result(msg)
        print_progress("Starting validation of reservoir case.")
        print_progress("reservoir_domain")
        validate_reservoir(reservoir_domain(case), reservoir_model(case), result)
        print_progress("state0")
        validate_dict(case.state0, model, result, :state0, :model)
        print_progress("parameters")
        validate_dict(case.parameters, model, result, :parameters, :model)
        print_progress("forces")
        validate_dt(case.dt, model, result)
        validate_forces(case.forces, case.dt, model, result)

        errors = oks = warnings = 0
        if info_level > 0
            print_result("Completed with $(length(result.messages)) messages.")
        end
        for (severity, prefix, msg) in result.messages
            if severity == OK
                oks += 1
                if info_level > 1
                    print_message(prefix, msg)
                end
            elseif severity == WARNING
                warnings += 1
                if info_level > 0
                    print_warning(prefix, msg)
                end
            elseif severity == ERROR
                errors += 1
                if info_level > -1
                    print_error(prefix, msg)
                end
            end
        end
        if info_level > 0
            if errors == 0 && warnings == 0
                print_success("Completed successfully.")
            else
                if errors > 0
                    if throw
                        error("Validation failed with $errors errors and $warnings warnings.")
                    end
                    print_error("Completed with $errors errors and $warnings warnings.")
                else
                    print_warning("Completed with no errors and $warnings warnings.")
                end
                println("Warnings indicate that values may be outside of expected ranges, but the case is still possible to simulate. Errors means that simulation will most likely fail.")
            end
        end
        return (errors == 0 && warnings == 0, result.messages)
    end

    function validate_dict(d, mm::MultiModel, result, name, model_name)
        for (k, m) in pairs(mm.models)
            haskey(d, k) || validation_error(result, name, "Subdict ismissing for submodel $model_name.")
            validate_dict(d[k], m, result, name, k)
        end
        return d
    end

    function validate_dict(d, model::SimulationModel, result, name, model_name)
        for (varname, varval) in pairs(d)
            vardef = get_variable(model, varname, throw = false)
            validate_variable(varval, vardef, model, result, varname, name, model_name)
        end
    end

    function validate_variable(varval, vardef::JutulVariables, model, result, varname, name, model_name)
        size(varval) == (Jutul.values_per_entity(model, vardef))
    end

    function validate_variable(varval, vardef::Nothing, model, result, varname, name, model_name)
        if varname == :ConnectionPressureDrop
            return
        end
        validation_message(result, name, "Variable $varname for $model_name was missing. May be unused or represent added variable.")
    end

    function validate_forces(forces::AbstractVector, dt, model, result)
        nforce = length(forces)
        ndt = length(dt)
        nforce == ndt || validation_error(result, "forces", "Number of forces ($nforce) and time steps do not match ($ndt).")
        for (i, force) in enumerate(forces)
            validate_forces(force, dt, model, result, i)
        end
    end

    function validate_forces(forces, dt, model, result, step = missing)
        for (k, v) in pairs(model.models)
            model_forces = get(forces, k, missing)
            ismissing(model_forces) && validation_warning(result, "forces", "Forces for submodel $k are missing at step $k.")
            if haskey(model_forces, :control)
                validate_facility_forces(model_forces, dt, v, result, step)
            end
        end
    end

    function validate_facility_forces(forces, dt, model, result, step)
        if ismissing(forces)
            return
        end
        ctrls = forces.control
        for (well, ctrl) in pairs(ctrls)
            if ctrl isa DisabledControl
                continue
            end
            target = ctrl.target
            info = JutulDarcy.well_target_information(target)
            if ismissing(step)
                msgname = "$well control"
            else
                msgname = "$well control at step $step"
            end
            if info.is_rate
                check(abs(target.value), :rate, result, msgname)
            elseif info.unit_type == :pressure
                check(target.value, :bhp, result, msgname)
            end
        end
    end

    function validate_dt(dt, model, result)
        if length(dt) == 0
            validation_warning(result, "dt", "Time step array is empty.")
        end
        num_warned = 0
        for (i, dt_i) in enumerate(dt)
            if dt_i <= 0.0
                num_warned += 1
                validation_error(result, "dt", "Time step at index $i is non-positive: $dt_i.")
            elseif dt_i < si_unit(:hour)
                num_warned += 1
                validation_warning(result, "dt", "Time step at index $i is very small: $dt_i seconds. Possible unit conversion error - time-steps should be given in seconds, not days.")
            end
            if num_warned > 10
                validation_warning(result, "dt", "More than 10 time steps are non-positive or very small. Further warnings will be suppressed.")
                break
            end
        end
    end

    function validate_reservoir(res::DataDomain, model, result)
        for k in keys(res)
            check(res[k], k, result, "reservoir_domain")
        end
    end

    # Utilities for checking
    function check(vals, label, result, name)
        ok = true
        messages = String[]
        fmt(x) = @sprintf("%.3g", x)
        if haskey(result.extrema, label)
            minval, maxval = extrema(vals)
            minval_def, maxval_def = result.extrema[label]
            if minval < minval_def 
                push!(messages, "$label has minimum value $(fmt(minval)), which is below the recommended minimum of $(fmt(minval_def)).")
                ok = false
            end
            if maxval > maxval_def
                push!(messages, "$label has maximum value $(fmt(maxval)), which is above the recommended maximum of $(fmt(maxval_def)).")
                ok = false
            end
        end
        if haskey(result.mean, label)
            meanval = mymean(vals)
            minmean, maxmean = result.mean[label]
            if meanval < minmean
                push!(messages, "$label has mean value $(fmt(meanval)), which is below the recommended minimum of $(fmt(minmean)).")
                ok = false
            end
            if meanval > maxmean
                push!(messages, "$label has mean value $(fmt(meanval)), which is above the recommended maximum of $(fmt(maxmean)).")
                ok = false
            end
        end

        nbad = 0
        firstbad = 0
        for (i, v) in enumerate(vals)
            if !isfinite(v)
                nbad += 1
                if firstbad == 0
                    firstbad = i
                end
            end
        end
        if nbad > 0
            push!(messages, "$label has $nbad non-finite or negative values. First occurrence at index $firstbad.")
            ok = false
        end

        if ok
            validation_message(result, name, "$label values passed validation.")
        else
            header = "$label - validation of values failed.\n"
            base_msg = join(messages, "\n")
            msg = header * base_msg
            u = get(result.default_units, label, missing)
            if !ismissing(u)
                msg = msg * "\nHint: Values are expected to be in unit $u. Check conversion."
            end
            if haskey(result.custom_messages, label)
                msg = msg *  "\nHint: " * result.custom_messages[label]
            end
            validation_warning(result, name, msg)
        end
    end

    function validation_message(result, prefix, msg)
        push!(result.messages, (OK, String(prefix), msg))
    end

    function validation_warning(result, prefix, msg)
        push!(result.messages, (WARNING, String(prefix), msg))
    end

    function validation_error(result, prefix, msg)
        if result.eager && result.throw
            error(msg)
        end
        push!(result.messages, (ERROR, String(prefix), msg))
    end

    function print_result(msg; color = :blue)
        print_result("Validation", msg, color = color)
    end

    function print_result(prefix, msg; color = :blue)
        jutul_message(prefix, msg, color = color)
    end

    function print_success(arg...)
        print_result(arg...; color = :green)
    end

    function print_message(arg...)
        print_result(arg...; color = :blue)
    end

    function print_warning(arg...)
        print_result(arg..., color = :yellow)
    end

    function print_error(arg...)
        print_result(arg..., color = :red)
    end

    function mymean(x)
        return sum(x) / length(x)
    end
end
