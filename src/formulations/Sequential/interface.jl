
function convert_to_sequential(model;
        avg_mobility = false,
        pressure = true,
        correction = missing,
        transport_scheme = :ppu
    )
    if pressure
        f = PressureFormulation()
        T_ctx = typeof(model.context)
        ctx = T_ctx(matrix_layout = EquationMajorLayout())
    else
        @assert ismissing(correction)
        f = TransportFormulation(schema = transport_scheme)
        ctx = model.context
        if model_or_domain_is_well(model)
            correction = :well
        else
            # correction = :volume
            correction = :density
        end
    end
    transport = !pressure
    seqmodel = SimulationModel(
        model.domain,
        model.system,
        data_domain = model.data_domain,
        formulation = f,
        context = ctx
    )
    for (skey, svar) in model.secondary_variables
        seqmodel.secondary_variables[skey] = svar
    end
    if pressure
        for (pkey, pvar) in model.primary_variables
            if pkey != :Pressure
                seqmodel.parameters[pkey] = pvar
            end
        end
        if avg_mobility
            mob = :PhaseMobilities
            if haskey(seqmodel.secondary_variables, mob)
                seqmodel.parameters[mob] = JutulDarcy.PhaseMobilities()
                delete!(seqmodel.secondary_variables, mob)
            end
        end
    end

    if transport
        vars = seqmodel.secondary_variables
        prm = seqmodel.parameters
        nph = number_of_phases(seqmodel.system)
        # has_mobility = haskey(vars, :PhaseMassMobilities)
        # has_mass_mobility = haskey(vars, :PhaseMassMobilities)
        #   # if correction == :volume
        #     if !has_mobility && !has_mass_mobility && !has_saturations
        #         correction = :density
        #     end
        # end has_saturations = haskey(vars, :Saturations)
        if correction == :density
            set_density_correction!(vars, prm, nph)
        elseif correction == :well
            set_well_correction!(vars, prm, nph)
        else
            correction == :volume || throw(ArgumentError("Unknown correction type $correction"))
            set_volume_correction!(vars, prm, nph)
        end
        push!(seqmodel.output_variables, :Pressure)
    end
    for (pkey, pvar) in model.parameters
        seqmodel.parameters[pkey] = pvar
    end
    for (pkey, pvar) in seqmodel.parameters
        if haskey(seqmodel.secondary_variables, pkey)
            delete!(seqmodel.parameters, pkey)
        end
    end
    for k in model.output_variables
        push!(seqmodel.output_variables, k)
    end
    if transport
        if model_or_domain_is_well(model)
            push!(seqmodel.output_variables, :PerforationTotalVolumetricFlux)
        else
            push!(seqmodel.output_variables, :TotalVolumetricFlux)
        end
    end
    unique!(seqmodel.output_variables)
    return seqmodel
end

function transport_set_corrected_phase_variable!(vars, k, nph)
    sym = Symbol(k, :Uncorrected)
    @assert haskey(vars, k) "Cannot correct $k, not present in variables: $(keys(vars))"
    vars[sym] = vars[k]
    vars[k] = TotalSaturationCorrectedVariable(sym, nph)
    return vars
end

function set_density_correction!(vars, prm, nph)
    if haskey(vars, :ShrinkageFactors)
        k = :ShrinkageFactors
    else
        k = :PhaseMassDensities
    end
    transport_set_corrected_phase_variable!(vars, k, nph)
    return vars
end

function set_volume_correction!(vars, prm, nph)
    has_mobility = haskey(vars, :PhaseMassMobilities)
    has_mass_mobility = haskey(vars, :PhaseMassMobilities)
    has_saturations = haskey(vars, :Saturations)


    if haskey(vars, :SurfaceMassMobilities)
        transport_set_corrected_phase_variable!(vars, :SurfaceMassMobilities, nph)
    end
    if has_mass_mobility
        k = :PhaseMassMobilities
    elseif has_saturations
        k = :Saturations
    end
    @assert haskey(vars, k) "Expected $k in $(keys(vars))"
    transport_set_corrected_phase_variable!(vars, k, nph)

    if k != :Saturations
        if haskey(prm, :StaticFluidVolume)
            pv_key = :StaticFluidVolume
        else
            pv_key = :FluidVolume
        end
        prm[:UncorrectedFluidVolume] = prm[pv_key]
        vars[pv_key] = TotalSaturationCorrectedScalarVariable(:UncorrectedFluidVolume)
    end
    return vars
end

function set_well_correction!(vars, prm, nph)
    if haskey(prm, :StaticFluidVolume)
        pv_key = :StaticFluidVolume
    else
        pv_key = :FluidVolume
    end
    prm[:UncorrectedFluidVolume] = prm[pv_key]
    vars[pv_key] = TotalSaturationCorrectedScalarVariable(:UncorrectedFluidVolume)
    return vars
end

function convert_to_sequential(model::MultiModel; pressure = true, kwarg...)
    ct = Vector{Jutul.CrossTermPair}()
    for ctp in model.cross_terms
        ctp = deepcopy(ctp)
        (; target, source, target_equation, source_equation, cross_term) = ctp
        target_is_cons = target_equation == :mass_conservation
        source_is_cons = source_equation == :mass_conservation
        if pressure && (target_is_cons || source_is_cons)
            if source_is_cons
                source_equation = :pressure
            end
            if target_is_cons
                target_equation = :pressure
            end
            if target_is_cons && source_is_cons && cross_term isa ReservoirFromWellFlowCT
                cross_term_reversed = PressureWellFromReservoirFlowCT(cross_term)
                ctp_reversed = Jutul.CrossTermPair(source, target, source_equation, target_equation, cross_term_reversed)
                push!(ct, ctp_reversed)
                cross_term = PressureReservoirFromWellFlowCT(cross_term)
            end
            ctp = Jutul.CrossTermPair(target, source, target_equation, source_equation, cross_term)

        end
        push!(ct, ctp)
    end
    smodel = convert_to_sequential(model[:Reservoir]; pressure = pressure, kwarg...)
    models = OrderedDict{Symbol, Any}()
    for (k, v) in pairs(model.models)
        if k == :Reservoir
            models[k] = smodel
        else
            if v.system isa MultiPhaseSystem
                v = convert_to_sequential(v; pressure = pressure, kwarg...)
            end
            models[k] = deepcopy(v)
        end
    end

    if isnothing(model.groups)
        g = nothing
    else
        g = copy(model.groups)
    end
    seqmodel = MultiModel(
        models,
        cross_terms = ct,
        groups = g,
        context = model.context,
        reduction = model.reduction
        )
    return seqmodel
end

function JutulDarcy.reservoir_linsolve(model::PressureModel, pname = :amg;
        solver = :bicgstab,
        rtol = 1e-3,
        kwarg...
    )
    if pname == :amg
        prec = default_psolve()
        lsolve = GenericKrylov(solver; preconditioner = prec, rtol = rtol, kwarg...)
    else
        lsolve = nothing
    end
    return lsolve
end
