function cpr_kernel_block_size(matrix)
    return max(matrix.minbatch, 128)
end

function synchronize_cpr_kernel(backend, event)
    if !isnothing(event)
        wait(event)
    end
    return nothing
end

@kernel function update_pressure_system_kernel!(
        pressure_values, @Const(rowptr), @Const(system_values),
        @Const(weights), number_of_components,
        represented_as_adjoint)
    row = @index(Global)
    @inbounds for position in rowptr[row]:(rowptr[row + 1] - 1)
        block = system_values[position]
        value = zero(eltype(pressure_values))
        for component in 1:number_of_components
            if represented_as_adjoint
                value += block[1, component]*weights[component, row]
            else
                value += block[component, 1]*weights[component, row]
            end
        end
        pressure_values[position] = value
    end
end

function update_pressure_system!(pressure_matrix::Jutul.StaticSparsityMatrixCSR,
        pressure_preconditioner,
        system_matrix::Jutul.StaticSparsityMatrixCSR,
        weights, context::Jutul.KernelAbstractionsContext, executor,
        ::Nothing)
    size(pressure_matrix) == size(system_matrix) || throw(DimensionMismatch(
        "CPR and reservoir matrices must have matching dimensions"))
    length(nonzeros(pressure_matrix)) == length(nonzeros(system_matrix)) ||
        throw(DimensionMismatch(
            "CPR and reservoir matrices must have matching sparsity"))
    backend = context.backend
    kernel! = update_pressure_system_kernel!(
        backend, cpr_kernel_block_size(pressure_matrix))
    event = kernel!(
        nonzeros(pressure_matrix), pressure_matrix.rowptr,
        nonzeros(system_matrix), weights, size(weights, 1),
        Jutul.represented_as_adjoint(Jutul.matrix_layout(context));
        ndrange = size(pressure_matrix, 1))
    synchronize_cpr_kernel(backend, event)
    return pressure_matrix
end

@inline function reduce_cpr_block(block, weights, entity,
        number_of_components, represented_as_adjoint)
    value = zero(eltype(weights))
    @inbounds for component in 1:number_of_components
        if represented_as_adjoint
            value += block[1, component]*weights[component, entity]
        else
            value += block[component, 1]*weights[component, entity]
        end
    end
    return value
end

@kernel function update_cprw_reservoir_kernel!(pressure_values,
        @Const(pressure_positions), @Const(rowptr),
        @Const(system_values), @Const(weights), number_of_components,
        represented_as_adjoint)
    row = @index(Global)
    @inbounds for position in rowptr[row]:(rowptr[row + 1] - 1)
        pressure_position = pressure_positions[position]
        block = system_values[position]
        pressure_values[pressure_position] = reduce_cpr_block(
            block, weights, row, number_of_components,
            represented_as_adjoint)
    end
end

@kernel function update_cprw_cross_kernel!(pressure_values,
        @Const(pressure_positions), @Const(system_positions),
        @Const(system_values), @Const(weight_entities),
        @Const(weights), number_of_components, represented_as_adjoint)
    index = @index(Global)
    @inbounds begin
        positions = system_positions[index]
        block = map(position -> system_values[position], positions)
        entity = weight_entities[index]
        pressure_values[pressure_positions[index]] = reduce_cpr_block(
            block, weights, entity, number_of_components,
            represented_as_adjoint)
    end
end

function update_pressure_system!(pressure_matrix::Jutul.StaticSparsityMatrixCSR,
        pressure_preconditioner,
        system_matrix::Jutul.StaticSparsityMatrixCSR,
        weights, context::Jutul.KernelAbstractionsContext, executor,
        well_reservoir_map)
    backend_data = well_reservoir_map.backend_data[]
    isnothing(backend_data) && throw(ArgumentError(
        "CPRW backend mappings were not initialized"))
    backend = context.backend
    block_size = cpr_kernel_block_size(pressure_matrix)
    number_of_components = size(weights, 1)
    represented_as_adjoint = Jutul.represented_as_adjoint(
        Jutul.matrix_layout(context))
    pressure_values = nonzeros(pressure_matrix)
    fill!(pressure_values, zero(eltype(pressure_values)))

    reservoir_kernel! = update_cprw_reservoir_kernel!(backend, block_size)
    event = reservoir_kernel!(pressure_values, backend_data.nzmap_11,
        system_matrix.rowptr, nonzeros(system_matrix), weights,
        number_of_components, represented_as_adjoint;
        ndrange = size(system_matrix, 1))
    synchronize_cpr_kernel(backend, event)

    cross_kernel! = update_cprw_cross_kernel!(backend, block_size)
    cross_terms = (
        (backend_data.nzmap_12, backend_data.map_12,
            well_reservoir_map.nzval_12, backend_data.cells_12),
        (backend_data.nzmap_21, backend_data.map_21,
            well_reservoir_map.nzval_21, backend_data.wells_21),
        (backend_data.nzmap_22, backend_data.map_22,
            well_reservoir_map.nzval_22, backend_data.wells_22)
    )
    for (pressure_positions, system_positions, system_values,
            weight_entities) in cross_terms
        count = length(pressure_positions)
        iszero(count) && continue
        event = cross_kernel!(pressure_values, pressure_positions,
            system_positions, system_values, weight_entities, weights,
            number_of_components, represented_as_adjoint;
            ndrange = count)
        synchronize_cpr_kernel(backend, event)
    end
    return pressure_matrix
end

@inline function cpr_weight_vector(accumulation, rhs, cell,
        pressure_scale, ::Val{N}, ::Val{unit_scaling}) where {N, unit_scaling}
    scalar_type = eltype(rhs)
    matrix = MMatrix{N, N, scalar_type}(undef)
    @inbounds for component in 1:N
        value = accumulation[component, cell]
        matrix[1, component] = value.partials[1]*pressure_scale
        for variable in 2:N
            matrix[variable, component] = value.partials[variable]
        end
    end
    weights = SMatrix{N, N, scalar_type}(matrix) \ SVector{N, scalar_type}(rhs)
    scale = if unit_scaling
        inv(norm(weights))
    else
        one(scalar_type)
    end
    return weights*scale
end

@kernel function true_impes_weights_kernel!(
        weights, @Const(accumulation), rhs,
        pressure_scale, number_of_cells, components::Val{N},
        unit_scaling) where N
    cell = @index(Global)
    if cell <= number_of_cells
        values = cpr_weight_vector(
            accumulation, rhs, cell, pressure_scale, components, unit_scaling)
        @inbounds for component in 1:N
            weights[component, cell] = values[component]
        end
    end
end

function update_true_impes_weights!(weights, accumulation, rhs,
        number_of_cells, number_of_components, pressure_scale, scaling,
        context::Jutul.KernelAbstractionsContext)
    backend = context.backend
    kernel! = true_impes_weights_kernel!(backend, 128)
    event = kernel!(weights, accumulation, rhs, pressure_scale,
        number_of_cells, Val(number_of_components), Val(scaling == :unit);
        ndrange = number_of_cells)
    synchronize_cpr_kernel(backend, event)
    return weights
end

@kernel function quasi_impes_weights_kernel!(
        weights, @Const(rowptr), @Const(colval),
        @Const(system_values), rhs, number_of_cells,
        offset, components::Val{N}, unit_scaling,
        represented_as_adjoint) where N
    cell = @index(Global)
    if cell <= number_of_cells
        row = cell + offset
        block = zero(eltype(system_values))
        @inbounds for position in rowptr[row]:(rowptr[row + 1] - 1)
            if colval[position] == row
                block = system_values[position]
                break
            end
        end
        matrix = represented_as_adjoint ? block : transpose(block)
        values = matrix \ SVector{N, eltype(rhs)}(rhs)
        scale = unit_scaling isa Val{true} ? inv(norm(values)) : one(eltype(rhs))
        @inbounds for component in 1:N
            weights[component, cell] = values[component]*scale
        end
    end
end

function update_quasi_impes_weights!(weights,
        system_matrix::Jutul.StaticSparsityMatrixCSR, rhs,
        number_of_cells, number_of_components, scaling, offset,
        context::Jutul.KernelAbstractionsContext)
    backend = context.backend
    kernel! = quasi_impes_weights_kernel!(backend, 128)
    event = kernel!(weights, system_matrix.rowptr, Jutul.colvals(system_matrix),
        nonzeros(system_matrix), rhs, number_of_cells, offset,
        Val(number_of_components), Val(scaling == :unit),
        Jutul.represented_as_adjoint(Jutul.matrix_layout(context));
        ndrange = number_of_cells)
    synchronize_cpr_kernel(backend, event)
    return weights
end

@kernel function analytical_cpr_weights_kernel!(
        weights, @Const(density), number_of_cells,
        number_of_components)
    cell = @index(Global)
    if cell <= number_of_cells
        @inbounds for component in 1:number_of_components
            weights[component, cell] = inv(Jutul.value(density[component, cell]))
        end
    end
end

function update_analytical_cpr_weights!(weights, model, state, rhs,
        number_of_cells, number_of_components, scaling,
        context::Jutul.KernelAbstractionsContext)
    map = Jutul.global_map(model.domain)
    density = Jutul.active_view(
        state.PhaseMassDensities, map; for_variables = false)
    kernel! = analytical_cpr_weights_kernel!(context.backend, 128)
    event = kernel!(weights, density, number_of_cells, number_of_components;
        ndrange = number_of_cells)
    synchronize_cpr_kernel(context.backend, event)
    return weights
end

@kernel function pressure_rhs_kernel!(
        pressure_rhs, @Const(full_rhs), @Const(weights), number_of_cells,
        number_of_components, block_size, forward_mode)
    cell = @index(Global)
    if cell <= length(pressure_rhs)
        value = zero(eltype(pressure_rhs))
        if cell <= number_of_cells
            if forward_mode
                @inbounds for component in 1:number_of_components
                    value += full_rhs[(cell - 1)*block_size + component]*
                        weights[component, cell]
                end
            else
                @inbounds value = full_rhs[(cell - 1)*block_size + 1]*
                    weights[1, cell]
            end
        end
        @inbounds pressure_rhs[cell] = value
    end
end

function update_p_rhs!(pressure_rhs, full_rhs, number_of_components,
        block_size, weights,
        pressure_matrix::Jutul.StaticSparsityMatrixCSR{
            <:Any, <:Any, <:Any, <:Any, <:Any,
            <:KernelAbstractions.Backend}, mode)
    backend = pressure_matrix.backend
    kernel! = pressure_rhs_kernel!(backend, 128)
    number_of_cells = length(full_rhs) ÷ block_size
    event = kernel!(pressure_rhs, full_rhs, weights, number_of_cells,
        number_of_components, block_size, mode == :forward;
        ndrange = length(pressure_rhs))
    synchronize_cpr_kernel(backend, event)
    return pressure_rhs
end
