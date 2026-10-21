# Included into JutulDarcy only by the isolated probe until full-run validation.
function cpu_convergence_supported(r, volume, properties, dt = 1.0)
    return r isa StridedMatrix{Float64} && volume isa StridedVector &&
        properties isa StridedMatrix && dt isa Float64 &&
        Threads.nthreads(:default) > 1 && length(volume) >= 16_384 &&
        typeof(Jutul.value(zero(eltype(volume)))) === Float64 &&
        typeof(Jutul.value(zero(eltype(properties)))) === Float64
end

function cpu_convergence_volume(volume)
    total = 0.0
    @batch minbatch=1000 reduction=((+, total),) for cell in eachindex(volume)
        @inbounds total += Jutul.value(volume[cell])
    end
    return total
end

function cpu_convergence_phase(r, volume, properties, dt, phase, ::Val{blackoil}) where blackoil
    maximum_residual = -Inf
    residual_sum = 0.0
    property_sum = 0.0
    @batch minbatch=1000 reduction=((max, maximum_residual), (+, residual_sum), (+, property_sum)) for cell in eachindex(volume)
        @inbounds begin
            residual = r[phase, cell]
            pore_volume = Jutul.value(volume[cell])
            property = Jutul.value(properties[phase, cell])
        end
        residual_sum += residual
        if blackoil
            property_sum += inv(property)
            maximum_residual = max(maximum_residual, abs(residual)/pore_volume)
        else
            property_sum += abs(property)
            maximum_residual = max(maximum_residual, dt*abs(residual)/(property*pore_volume))
        end
    end
    return maximum_residual, residual_sum, property_sum/length(volume)
end

function cpu_convergence_errors(r, volume, properties, dt, ::Val{N}, blackoil::Val{BO}, reference_density=nothing) where {N, BO}
    size(r, 2) == size(properties, 2) == length(volume) || throw(DimensionMismatch())
    size(r, 1) >= N && size(properties, 1) >= N || throw(DimensionMismatch())
    total_volume = cpu_convergence_volume(volume)
    metrics = ntuple(phase -> cpu_convergence_phase(r, volume, properties, dt, phase, blackoil), Val(N))
    cnv = ntuple(Val(N)) do phase
        maximum_residual, _, property_mean = metrics[phase]
        BO ? property_mean*dt*maximum_residual/reference_density[phase] : maximum_residual
    end
    mb = ntuple(Val(N)) do phase
        _, residual_sum, property_mean = metrics[phase]
        BO ? property_mean*dt*abs(residual_sum)/(reference_density[phase]*total_volume) :
            (dt/total_volume)*abs(residual_sum)/property_mean
    end
    return cnv, mb
end
