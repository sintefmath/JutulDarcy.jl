struct EquilibriumRegion{R}
    datum_pressure::R
    datum_depth::R
    woc::R
    goc::R
    wgc::R
    pc_woc::R
    pc_goc::R
    pc_wgc::R
    density_function::Union{Function, Missing}
    composition_vs_depth::Union{Function, Missing}
    rs_vs_depth::Union{Function, Missing}
    rv_vs_depth::Union{Function, Missing}
    temperature_vs_depth::Union{Function, Missing}
    cells::Union{Vector{Int}, Missing}
    pvtnum::Int
    satnum::Int
    kwarg::Any
end

"""
    EquilibriumRegion(model, p_datum = missing,
        datum_depth = missing;
        woc = NaN,
        goc = NaN,
        wgc = NaN,
        pc_woc = 0.0,
        pc_goc = 0.0,
        pc_wgc = 0.0,
        temperature = missing,
        temperature_vs_depth = missing,
        composition = missing,
        composition_vs_depth = missing,
        liquid_composition = missing,
        liquid_composition_vs_depth = missing,
        vapor_composition = missing,
        vapor_composition_vs_depth = missing,
        density_function = missing,
        rs = 0.0,
        rs_vs_depth = missing,
        rv = 0.0,
        rv_vs_depth = missing,
        pvtnum = 1,
        satnum = 1,
        cells = missing,
        kwarg...
    )

Set uip equilibriation region for a reservoir model. The region is defined by
the datum pressure and depth, and the water-oil, gas-oil, and water-gas
contacts. The contacts will be used to determine the phase distribution and
initial pressure in the region. The region can be further specified by
temperature, composition, density, and rs/rv functions. Most entries can either
be specified as a function of depth or as a constant value. Additional keyword
arguments are passed onto the `equilibriate_state` function.
"""
function EquilibriumRegion(model::Union{SimulationModel, MultiModel}, p_datum = missing,
        datum_depth = missing;
        woc = NaN,
        goc = NaN,
        wgc = NaN,
        pc_woc = 0.0,
        pc_goc = 0.0,
        pc_wgc = 0.0,
        temperature = missing,
        temperature_vs_depth = missing,
        composition = missing,
        composition_vs_depth = missing,
        liquid_composition = missing,
        liquid_composition_vs_depth = missing,
        vapor_composition = missing,
        vapor_composition_vs_depth = missing,
        density_function = missing,
        rs = 0.0,
        rs_vs_depth = missing,
        rv = 0.0,
        rv_vs_depth = missing,
        pvtnum = 1,
        satnum = 1,
        cells = missing,
        kwarg...
    )
    function handle_depth_or_value(x, xvd, xname; num = 0, is_density = false)
        if ismissing(x) && ismissing(xvd)
            error("Either $xname or $(xname)_vs_depth must be provided")
        elseif !ismissing(xvd)
            # Note: vs depth tables take precedence
            xvd isa Function || error("$(xname)_vs_depth must be a function on the form z -> val")
            F = xvd
        else
            if num > 0
                length(x) == num || error("Length of $xname must be $num")
            else
                x isa Real || error("$xname must be a real number")
            end
            F = z -> x
        end
        return F
    end
    model = reservoir_model(model)
    reservoir = reservoir_domain(model)
    rmesh = physical_representation(reservoir)
    # Checks for various compositions vs depth
    sys = model.system
    is_compositional = sys isa MultiPhaseCompositionalSystemLV
    if ismissing(density_function)
        if has_disgas(sys)
            rs_vs_depth = handle_depth_or_value(rs, rs_vs_depth, "rs")
        end
        if has_vapoil(sys)
            rv_vs_depth = handle_depth_or_value(rv, rv_vs_depth, "rv")
        end
        if is_compositional
            n_hc = MultiComponentFlash.number_of_components(sys.equation_of_state)
            if ismissing(composition_vs_depth) && ismissing(composition)
                x_vs_depth = handle_depth_or_value(liquid_composition, liquid_composition_vs_depth, "liquid_composition", num = n_hc)
                y_vs_depth = handle_depth_or_value(vapor_composition, vapor_composition_vs_depth, "vapor_composition", num = n_hc)
                composition_vs_depth = z -> ifelse(z > goc, x_vs_depth(z), y_vs_depth(z))
            else
                composition_vs_depth = handle_depth_or_value(composition, composition_vs_depth, "composition", num = n_hc)
            end
        end
    end
    vars = Jutul.get_variables_by_type(model, :all)
    if haskey(vars, :Temperature)
        temperature_vs_depth = handle_depth_or_value(temperature, temperature_vs_depth, "temperature")
    end
    # Validate the cells
    nc = number_of_cells(rmesh)
    if ismissing(cells)
        cells = collect(1:nc)
    else
        cells = collect(cells)
        for c in cells
            if c < 1 || c > nc
                error("Cell index $c out of range 1:$nc")
            end
        end
    end
    return EquilibriumRegion(
        p_datum,
        datum_depth,
        woc,
        goc,
        wgc,
        pc_woc,
        pc_goc,
        pc_wgc,
        density_function,
        composition_vs_depth,
        rs_vs_depth,
        rv_vs_depth,
        temperature_vs_depth,
        cells,
        pvtnum,
        satnum,
        kwarg
    )
end

function Base.show(io::IO, eql::EquilibriumRegion{R}) where R
    nc = length(eql.cells)
    println(io, "EquilibriumRegion{$R} for $(nc) cells with datum pressure $(eql.datum_pressure) Pa at depth $(eql.datum_depth) m")
    print(io, "Water-Oil Contact: $(eql.woc) m, Gas-Oil Contact: $(eql.goc) m, Water-Gas Contact: $(eql.wgc) m\n")
end

"""
    setup_reservoir_state(model, <keyword arguments>)
    # Ex: For immiscible two-phase
    setup_reservoir_state(model, Pressure = 1e5, Saturations = [0.2, 0.8])

Convenience constructor that initializes a state for a `MultiModel` set up using
[`setup_reservoir_model`](@ref). The main convenience over [`setup_state`](@ref)
is only the reservoir initialization values need be provided: wells are
automatically initialized from the connected reservoir cells.

As an alternative to passing keyword arguments, a `Dict{Symbol, Any}` instance
can be sent in as a second, non-keyword argument.
"""
function setup_reservoir_state(model::MultiModel, equil::Union{Missing, Vector, EquilibriumRegion} = missing; kwarg...)
    rmodel = reservoir_model(model)
    pvars = collect(keys(Jutul.get_primary_variables(rmodel)))
    res_state = setup_reservoir_state(rmodel, equil; kwarg...)
    # Next, we initialize the wells.
    init = Dict(:Reservoir => res_state)
    perf_subset(v::AbstractVector, i) = v[i]
    perf_subset(v::AbstractMatrix, i) = v[:, i]
    perf_subset(v, i) = v
    is_thermal = model_is_thermal(rmodel)
    for k in keys(model.models)
        if k == :Reservoir
            # Already done
            continue
        end
        W = model.models[k]
        if W.domain isa WellGroup
            # Facility or well group
            if !is_thermal
                init_w = setup_state(W; TotalSurfaceMassRate = 0.0)
            else
                T0 = convert_to_si(20.0, :Celsius)
                init_w = setup_state(W; TotalSurfaceMassRate = 0.0, SurfaceTemperature = T0)
            end
        else
            # Wells
            init_w = Dict{Symbol, Any}()
            W = model.models[k]
            wg = physical_representation(W.domain)
            res_c = wg.perforations.reservoir
            if wg isa MultiSegmentWell
                init_w[:TotalMassFlux] = 0.0
            end
            c = map_well_nodes_to_reservoir_cells(wg, rmodel.data_domain)
            for pk in pvars
                pv = res_state[pk]
                init_w[pk] = perf_subset(pv, c)
            end
        end
        init[k] = init_w
    end
    state = setup_state(model, init)
    return state
end

function setup_reservoir_state(model, init::AbstractDict; kwarg...)
    if haskey(init, :Reservoir) && model isa MultiModel
        # Could be output from a previous call to the same routine
        init = init[:Reservoir]
    end
    return setup_reservoir_state(model; pairs(init)..., kwarg...)
end

function setup_reservoir_state(
        rmodel::SimulationModel,
        equil_regs::Union{Missing, Vector, EquilibriumRegion} = missing;
        kwarg...
    )
    if ismissing(equil_regs)
        init = kwarg
    else
        if equil_regs isa EquilibriumRegion
            equil_regs = [equil_regs]
        end
        inits = map(equil -> equilibriate_state(rmodel, equil), equil_regs)
        if length(inits) == 1
            init = only(inits)
        else
            # Handle multiple regions by merging each init
            inits_cells = map(x -> x.cells, equil_regs)
            init = Dict{Symbol, Any}()
            nc = number_of_cells(rmodel.domain)
            touched = [false for i in 1:nc]
            for (k, v) in first(inits)
                if v isa AbstractVector
                    init[k] = zeros(nc)
                else
                    @assert v isa AbstractMatrix
                    init[k] = zeros(size(v, 1), nc)
                end
            end
            for (subinit, cells) in zip(inits, inits_cells)
                for c in cells
                    if touched[c]
                        @warn "Equils overlap for cell $c?"
                    end
                    touched[c] = true
                end
                for (k, v) in subinit
                    fill_subinit!(init[k], cells, v)
                end
            end
            @assert all(touched) "Some cells are not initialized by equil: $(findall(!, touched))"
        end
        for (k, v) in kwarg
            init[k] = v
        end
    end

    pvars = collect(keys(Jutul.get_primary_variables(rmodel)))
    svars = collect(keys(Jutul.get_secondary_variables(rmodel)))
    np = length(pvars)
    found = Symbol[]
    res_init = Dict{Symbol, Any}()
    for (k, v) in init
        I = findfirst(isequal(k), pvars)
        if eltype(v)<:AbstractFloat
            if !all(isfinite, v)
                jutul_message("setup_reservoir_state", "Non-finite entries found in initializer for $k.", color = :red)
            end
        end
        if isnothing(I)
            if !(k in svars)
                jutul_message("setup_reservoir_state", "Received primary variable $k, but this is not known to reservoir model.")
            end
        else
            push!(found, k)
        end
        res_init[k] = v
    end
    tc = get(rmodel.primary_variables, :TracerConcentrations, nothing)
    if  !isnothing(tc) && !haskey(res_init, :TracerConcentrations)
        # Tracers are usually safe to default = 0
        res_init[:TracerConcentrations] = Jutul.default_values(rmodel, tc)
        push!(found, :TracerConcentrations)
    end
    handle_alternate_primary_variable_spec!(res_init, found, rmodel, rmodel.system)
    if length(found) != length(pvars)
        missing_primary_variables = setdiff(pvars, found)
        @warn "Not all primary variables were initialized for reservoir model." missing_primary_variables
    end
    return setup_state(rmodel, res_init)
end

function handle_alternate_primary_variable_spec!(res_init, found, rmodel, system)
    # Internal utility to handle non-trivial specification of primary variables
    return res_init
end
