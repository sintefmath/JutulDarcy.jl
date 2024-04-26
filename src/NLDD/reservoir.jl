function boundary_discretization(storage_g, model_g, submodel::MultiModel, storage)
    tag = Jutul.submodel_ad_tag(submodel, :Reservoir)
    return boundary_discretization(storage_g, model_g.models.Reservoir, submodel.models.Reservoir, storage.Reservoir, tag)
end

function boundary_discretization(storage_g, model_g, submodel, storage, tag = nothing)
    d = submodel.domain
    M = global_map(d)
    mass_flow = d.discretizations.mass_flow
    mf = reverse_bnd(mass_flow, M)

    # Deal with storage for the reversed fluxes
    nf = mf.conn_pos[end]-1
    mc = storage.equations.mass_conservation
    base_flux = Jutul.get_entries(mc.half_face_flux_cells)
    nph = size(base_flux, 1)
    T = eltype(base_flux)
    flux = zeros(T, nph, nf)
    # tagof(v::ForwardDiff.Dual{T}) where T = T
    # flux = JutulAutoDiffCache()
    eq = submodel.equations[:mass_conservation]
    ne = Jutul.number_of_equations_per_entity(submodel, eq)
    flux = CompactAutoDiffCache(ne, nf, submodel, entity = Cells(), tag = tag)

    lsys = storage_g.LinearizedSystem[1, 1]
    jac = lsys.jac
    context = submodel.context
    @assert matrix_layout(context) == matrix_layout(model_g.context)
    nu, ne, np = Jutul.ad_dims(flux)
    nu = number_of_cells(model_g.domain)

    cells = zeros(Int64, 0)
    for f_ix in 1:nf
        cd_f = mf.conn_data[f_ix]
        # f = cd_f.face
        # Cell outside active domain - corresponds to col ix
        self = cd_f.self
        @assert Jutul.cell_is_boundary(self, M)
        # Cell inside active domain - corresponds to row ix
        other = cd_f.other
        # Note: Might be wrong.
        # target = Jutul.global_cell(self, M)
        #source = Jutul.global_cell(other, M)
        source = Jutul.global_cell(self, M)
        target = Jutul.global_cell(other, M)
        # target = row
        # source = column
        push!(cells, self)
        # @info "$target -> $source"
        for e in 1:ne
            for d = 1:np
                # pos = Jutul.find_jac_position(jac, other_i, active_cell_i, e, d, nu, nu, ne, np, context)
                pos = Jutul.find_jac_position(jac, target, source, 0, 0, 0, 0, e, d, nu, nu, ne, np, context)
                Jutul.set_jacobian_pos!(flux, f_ix, e, d, pos)
            end
        end
    end
    @assert nf == length(mf.conn_data)
    return (flux = flux, mass_flow = mf, cells = cells)
end

update_aspen_coupling!(storage, substorage, submodel::MultiModel, state, p_i) = update_aspen_coupling!(storage, substorage.Reservoir, submodel.models.Reservoir, state.Reservoir, p_i)

function update_aspen_coupling!(storage, substorage, submodel, state, p_i)
    outer_state = substorage.state
    bnd = storage.boundary_discretizations[p_i]

    flow_disc = bnd.mass_flow
    flux = bnd.flux
    cells = bnd.cells

    state = handle_subdomain_boundary_variable_switching!(state, substorage, submodel, cells)
    eq = submodel.equations[:mass_conservation]
    # Put flux into right form and localize AD
    flux_c = Jutul.get_entries(flux)
    N = size(flux_c, 1)
    T = eltype(flux_c)
    flux_static = reinterpret(SVector{N, T}, flux_c)
    dt = NaN
    state_c = local_ad(state, 1, T)
    aspen_update_half_face_flux_tpfa!(flux_static, eq, state_c, submodel, dt, flow_disc, Cells())
end

function check_aspen_coupling!(storage, model, bnd, block_no)
    if model isa MultiModel
        rmodel = model[:Reservoir]
        rstorage = storage.Reservoir
    else
        rmodel = model
        rstorage = storage
    end
    disc = rmodel.domain.discretizations.mass_flow
    inner = rstorage.equations.mass_conservation

    conn_pos_outer = disc.conn_pos
    conn_data_outer = disc.conn_data
    flux_outer = inner.half_face_flux_cells.entries

    flow_disc = bnd.mass_flow
    flux = bnd.flux.entries
    cells = bnd.cells
    conn_data = flow_disc.conn_data
    conn_pos = flow_disc.conn_pos

    Nc = length(conn_pos)-1
    # @info "" size(flux_outer) size(conn_data_outer) conn_pos_outer flux_outer'
    # @assert size(flux_outer, 2) == conn_pos_outer[end] "$(size(flux_outer)) $(conn_pos_outer[end])"
    for c in 1:Nc
        for hf in conn_pos[c]:(conn_pos[c+1]-1)
            f = conn_data[hf].face
            q_bnd = flux[:, hf]

            outer_pos = 0
            for c_i in eachindex(conn_data_outer)
                if conn_data_outer[c_i].face == f
                    outer_pos = c_i 
                    break
                end
            end
            @assert outer_pos > 0
            c_o = conn_data_outer[outer_pos]
            @assert c_o.self == conn_data[hf].other
            @assert c_o.other == conn_data[hf].self
            q = flux_outer[:, outer_pos]
            for (in, out, i) in zip(value(q), -value(q_bnd), 1:length(q))
                err = abs(in - out)/abs(out)
                if err > 1e-12
                    msg = "|e|=$err: Block $block_no #$c/$(Nc) for face $f = $in != $out"
                    @info msg
                    if i == length(q)
                        @info "With partials" q q_bnd
                    end
                    # error()
                end
            end
        end
    end
    # error()
end

function aspen_update_half_face_flux_tpfa!(hf_cells::Union{AbstractArray{SVector{N, T}}, AbstractVector{T}}, eq, state::S, model, dt, flow_disc, ::Cells) where {T, N, S<:Jutul.LocalStateAD}
    # This is way too coupled to internal Jutul guts and should be replaced by something more general.
    conn_data = flow_disc.conn_data
    conn_pos = flow_disc.conn_pos
    Te = eltype(hf_cells)
    nhf = conn_pos[end]-1
    @tic "flux (cells)" for i in 1:nhf
        @inbounds (; self, other, face, face_sign) = conn_data[i]
        state_c = Jutul.new_entity_index(state, self)
        @inbounds hf_cells[i] = Jutul.face_flux!(zero(Te), self, other, face, face_sign, eq, state_c, model, dt, flow_disc)
    end
end


function reverse_bnd(disc::TwoPointPotentialFlowHardCoded, map::Jutul.FiniteVolumeGlobalMap)
    has_grav = disc.gravity

    conn_pos = disc.conn_pos
    conn_data = disc.conn_data

    is_bnd = map.cell_is_boundary

    nc = length(conn_pos)-1
    cell_is_near_bnd = [false for _ in 1:nc]
    for c in conn_data
        if is_bnd[c.other]
            c_local = Jutul.interior_cell(c.self, map)
            cell_is_near_bnd[c_local] = true
        end
    end
    cells = findall(cell_is_near_bnd)
    nc_new = length(cells)
    conn_data_new = similar(conn_data, 0)
    conn_pos_new = zeros(Int64, nc_new+1)
    ctr = 1
    conn_pos_new[1] = ctr
    for (i, c_i) in enumerate(cells)
        for f in conn_pos[c_i]:conn_pos[c_i+1]-1
            cd = conn_data[f]
            if is_bnd[cd.other]
                push!(conn_data_new, reverse_connection(cd))
                # @info "renumb" cd.other reverse_connection(cd)
                ctr += 1
            end
        end
        # @info "$i" ctr - conn_pos_new[i]
        conn_pos_new[i + 1] = ctr
    end
    # @info "All done" conn_pos_new conn_data_new ctr
    tp_flow = TwoPointPotentialFlowHardCoded{typeof(conn_pos_new), typeof(conn_data_new)}(has_grav, conn_pos_new, conn_data_new)
    return tp_flow
end

function reverse_connection(conn::T) where T
    D = OrderedDict()
    for k in keys(conn)
        if k == :self
            newval = conn.other
        elseif k == :other
            newval = conn.self
        elseif k == :face_sign
            newval = -conn.face_sign
        else
            newval = conn[k]
        end
        D[k] = newval
    end
    return convert_to_immutable_storage(D)::T
end

function black_oil_primary_buffers(storage, model, sys::JutulDarcy.BlackOilVariableSwitchingSystem)
    nc = Jutul.count_active_entities(model.domain, Cells())
    # rs, rv, sg
    return zeros(3, nc)
end

function black_oil_primary_buffers(storage, model, sys)
    return nothing
end


function black_oil_primary_buffers(storage, model)
    return black_oil_primary_buffers(storage, model, model.system)
end

function black_oil_primary_buffers(storage, model::MultiModel)
    buf = Dict{Symbol, Any}()
    groups = model.groups
    if isnothing(groups)
        offset = 0
        ix = 1
        for (k, m) in pairs(model.models)
            buf[k] = black_oil_primary_buffers(storage[k], m)
        end
    else
        ug = unique(groups)
        buffers = Vector{Matrix{Float64}}(undef, length(ug))
        for g in ug
            offset = 0
            ix = 1
            for (k, m) in pairs(model.models)
                if groups[ix] == g
                    buf[k] = black_oil_primary_buffers(storage[k], m, m.system)
                end
                ix += 1
            end
        end
    end
    return buf
end

function handle_subdomain_boundary_variable_switching!(state, substorage, submodel, cells)
    return state
end


import JutulDarcy: s_removed
import Jutul: replace_value
function handle_subdomain_boundary_variable_switching!(state, substorage, model::JutulDarcy.StandardBlackOilModel, cells)
    outer_state = substorage.state
    m = model.domain.global_map
    sys = model.system
    rs_max = sys.rs_max
    rv_max = sys.rv_max
    disgas = JutulDarcy.has_disgas(sys)
    vapoil = JutulDarcy.has_vapoil(sys)
    if haskey(model.secondary_variables, :Rs)
        boreg = model[:Rs].regions
    elseif haskey(model.secondary_variables, :Rv)
        boreg = model[:Rv].regions
    else
        boreg = nothing
    end
    for c in cells
        gc = Jutul.global_cell(c, m)
        reg_i = JutulDarcy.region(boreg, c)
        phases = state.PhaseState[c]
        g_phases = outer_state.PhaseState[gc]
        if phases != g_phases
            x = state.BlackOilUnknown[c]
            Sw = state.ImmiscibleSaturation[c]
            sg = value(state.Saturations[3, c])
            @inbounds if g_phases == JutulDarcy.OilAndGas
                x_val = sg
            else
                if g_phases == JutulDarcy.GasOnly
                    x_val = value(state.Rv[c])
                    Sg = 1 - Sw
                    So = 0
                else
                    x_val = value(state.Rs[c])
                    Sg = 0
                    So = 1 - Sw
                end
                # Also need to handle saturations - set to values but no
                # derivatives with respect to X
                state.Saturations[2, c] = So
                state.Saturations[3, c] = Sg
            end
            val = replace_value(x.val, x_val)
            state.BlackOilUnknown[c] = BlackOilX(val, g_phases)
            if g_phases == JutulDarcy.OilAndGas
                p = state.Pressure[c]
                # p = outer_state.Pressure[gc]
                @inbounds if disgas
                    rstab = JutulDarcy.table_by_region(model.system.rs_max, reg_i)
                    rs = rstab(p)
                    rs0 = value(state.Rs[c])
                    # Rs with value of old value and derivatives of new one
                    new_rs = replace_value(rs, rs0)
                    # new_rs = min(rs0/value(rs), 1.0)*rs
                    state.Rs[c] = new_rs
                end
                @inbounds if vapoil
                    rvtab = table_by_region(model.system.rs_max, reg_i)
                    rv = rvtab(p)
                    rv0 = value(state.Rv[c])
                    # new_rv = min(rv0/value(rv), 1.0)*rv
                    new_rv = replace_value(rv, rv0)
                    state.Rv[c] = new_rv
                end
            elseif g_phases == JutulDarcy.OilOnly
                # Then Rs as free variable should be ok
                @inbounds state.Rs[c] = val
                if vapoil
                    @inbounds state.Rv[c] = 0.0
                end
            else
                @inbounds state.Rv[c] = val
                if disgas
                    @inbounds state.Rs[c] = 0.0
                end
            end
        end
    end
    reupdate_secondary_variables_bo_specific!(state, model, substorage.variable_definitions.secondary_variables, cells)
    return state
end

function reupdate_secondary_variables_bo_specific!(state, model, secondaries, cells)
    skipped = (:Rs, :Rv, :Saturations, :TotalMasses, :FluidVolume)
    for (s, var) in pairs(secondaries)
        if s in skipped
            continue
        end
        v = state[s]
        update_secondary_variable!(v, var, model, state, cells)
    end
    return state
end

function reservoir_change_buffers(storage, model)
    n = number_of_cells(model.domain)
    nph = number_of_phases(model.system)
    buf = Dict(
        :Saturations => zeros(nph, n),
        :PhaseMobilities => zeros(nph, n),
        :Pressure => zeros(n)
        )
    if model isa JutulDarcy.CompositionalModel
        ncomp = JutulDarcy.number_of_components(model.system) - Int(JutulDarcy.has_other_phase(model.system))
        buf[:OverallMoleFractions] = zeros(ncomp, n)
    end
    buf = NamedTuple(pairs(buf))
    return buf::NamedTuple
end

function reservoir_change_buffers(storage, model::MultiModel)
    return reservoir_change_buffers(storage[:Reservoir], model[:Reservoir])
end

function store_reservoir_change_buffer!(buf, sim, cfg)
    tol_s = cfg[:solve_tol_saturations]
    tol_p = cfg[:solve_tol_pressure]
    tol_λ = cfg[:solve_tol_mobility]
    tol_z = cfg[:solve_tol_composition]
    model = sim.model
    state = sim.storage.state
    if model isa MultiModel
        model = model[:Reservoir]
        state = state.Reservoir
    end

    if !isnothing(tol_s) || !isnothing(tol_p) || !isnothing(tol_λ) || !isnothing(tol_z)
        store_reservoir_change_buffer_inner!(buf, model, state)
    end
end

function get_nldd_solution_change_tolerances(cfg)
    tol_s = cfg[:solve_tol_saturations]
    tol_p = cfg[:solve_tol_pressure]
    tol_mob = cfg[:solve_tol_mobility]
    tol_z = cfg[:solve_tol_composition]

    has_s = !isnothing(tol_s)
    has_p = !isnothing(tol_p)
    has_mob = !isnothing(tol_mob)
    has_z = !isnothing(tol_z)
    if has_s || has_p || has_mob || has_z
        out = (
            s = tol_s,
            p = tol_p,
            mob = tol_mob,
            z = tol_z
        )
    else
        out = nothing
    end
    return out
end

function check_if_subdomain_needs_solving(buf, sim, cfg, iteration)
    tol = get_nldd_solution_change_tolerances(cfg)
    model = sim.model
    has_wells = model isa MultiModel && length(model.models) > 1
    if cfg[:always_solve_wells] && has_wells
        return true
    end
    if isnothing(tol)
        do_solve = true
    else
        if iteration == 1
            do_solve = !cfg[:solve_tol_first_newton]
        else
            state = sim.storage.state
            do_solve = check_inner(buf, model, state, tol)
        end
    end
    return do_solve
end

function check_inner(buf, model, state, tol)
    do_solve = false
    do_solve = do_solve || check_subdomain_change_inner(buf, model, state, :Saturations, tol.s, :abs)
    do_solve = do_solve || check_subdomain_change_inner(buf, model, state, :Pressure, tol.p, :abs)
    do_solve = do_solve || check_subdomain_change_inner(buf, model, state, :PhaseMobilities, tol.mob, :relsum)
    do_solve = do_solve || check_subdomain_change_inner(buf, model, state, :OverallMoleFractions, tol.z, :abs)
    return do_solve
end

function store_reservoir_change_buffer_inner!(buf, model, state)
    function buf_transfer!(out, in)
        for i in eachindex(out)
            @inbounds out[i] = value(in[i])
        end
    end
    for (k, v) in pairs(buf)
        vals = state[k]
        buf_transfer!(v, vals)
    end
end

function check_subdomain_change_inner(buf, model::MultiModel, state, f, tol, sum_t)
    return check_subdomain_change_inner(buf, model[:Reservoir], state[:Reservoir], f, tol, sum_t)
end

function check_subdomain_change_inner(buf, model::SimulationModel, state, f, tol, sum_t)
    if isnothing(tol)
        return false # No tolerance - no need to check
    else
        current = state[f]
        old = buf[f]
        if sum_t == :relsum
            current::AbstractMatrix
            return subdomain_delta_relsum(current, old, tol)
        else
            @assert sum_t == :abs
            return subdomain_delta_absolute(current, old, tol)
        end
    end
end

function subdomain_delta_absolute(current, old, tol)
    for i in eachindex(current)
        c = current[i]
        o = old[i]
        d = abs(value(c) - o)
        if d > tol
            return true
        end
    end
    return false
end

function subdomain_delta_relsum(current, old, tol)
    n = size(current, 1)
    for cell in axes(current, 2)
        mob_total = 0.0
        for i in axes(old, 1)
            mob_total += old[i, cell]/n
        end
        for i in axes(old, 1)
            c = current[i, cell]
            o = old[i, cell]
            d = abs(value(c) - o)/mob_total
            if d > tol
                return true
            end
        end
    end
    return false
end
