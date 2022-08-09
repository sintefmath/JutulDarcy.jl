export DeckViscosity, DeckShrinkage

function check_regions(regions, t)
    if !isnothing(regions)
        regions<:AbstractVector
        @assert maximum(regions) <= length(t)
        @assert minimum(regions) > 0
    end
end


region(pv::DeckPhaseVariables, cell) = region(pv.regions, cell)
region(r::AbstractVector, cell) = @inbounds r[cell]
region(::Nothing, cell) = 1

function tab_by_region(pvt, reg)
    return pvt.tab[reg]
end


@jutul_secondary function update_as_secondary!(mu, μ::DeckViscosity, model, param, Pressure)
    pvt, reg = μ.pvt, μ.regions
    if false
        @tullio mu[ph, i] = viscosity(pvt[ph], reg, Pressure[i], i)
    else
        tb = minbatch(model.context)
        nph, nc = size(mu)
        for ph in 1:nph
            pvt_ph = pvt[ph]
            @batch minbatch = tb for i in 1:nc
                p = Pressure[i]
                @inbounds mu[ph, i] = viscosity(pvt_ph, reg, p, i)
            end
        end
    end
end

@jutul_secondary function update_as_secondary!(rho, ρ::DeckDensity, model, param, Pressure)
    rhos = reference_densities(model.system)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    if false
        function do_stuff!(rho, rhos, reg, Pressure)
            @tullio rho[ph, i] = rhos[ph]*shrinkage(pvt[ph], reg, Pressure[i], i)
        end
        do_stuff!(rho, rhos, reg, Pressure)
    else
        tb = minbatch(model.context)
        nph, nc = size(rho)
        for ph in 1:nph
            rhos_ph = rhos[ph]
            pvt_ph = pvt[ph]
            @batch minbatch = tb for i in 1:nc
                p = Pressure[i]
                @inbounds rho[ph, i] = rhos_ph*shrinkage(pvt_ph, reg, p, i)
            end
        end
    end
    # @tullio rho[ph, i] = rhos[ph]*shrinkage(pvt[ph], reg, Pressure[i], i)
end


@jutul_secondary function update_as_secondary!(b, ρ::DeckShrinkageFactors, model, param, Pressure)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    tb = minbatch(model.context)
    nph, nc = size(b)
    for ph in 1:nph
        pvt_ph = pvt[ph]
        @batch minbatch = tb for i in 1:nc
            p = Pressure[i]
            @inbounds b[ph, i] = shrinkage(pvt_ph, reg, p, i)
        end
    end
end

@jutul_secondary function update_as_secondary!(pv, Φ::LinearlyCompressiblePoreVolume, model, param, Pressure)
    vol = Φ.volume
    c_r = Φ.expansion
    p_r = Φ.reference_pressure
    for i in eachindex(pv)
        p = Pressure[i]
        x = c_r*(p-p_r)
        mult = 1 + x + 0.5*(x^2)
        pv[i] = vol[i]*mult
    end
end

# struct DeckRelativePermeability <: DeckPhaseVariables
#     sat::Tuple
#     regions
#     function DeckRelativePermeability(sat; regions = nothing)
#         check_regions(regions, sat)
#         new(Tuple(sat), regions)
#     end
# end

# @jutul_secondary function update_as_secondary!(kr, kr_def::DeckRelativePermeability, model, param, Saturations)
#     @tullio kr[ph, i] = relative_permeability(get_sat(kr_def, ph, i), Saturations[ph, i])
# end

# struct DeckCapillaryPressure <: DeckPhaseVariables
#     sat::Tuple
#     regions
#     function DeckRelativePermeability(sat; regions = nothing)
#         check_regions(regions, sat)
#         new(Tuple(sat), regions)
#     end
# end

# degrees_of_freedom_per_entity(model, sf::DeckCapillaryPressure) = number_of_phases(model.system) - 1


# @jutul_secondary function update_as_secondary!(kr, kr_def::DeckCapillaryPressure, model, param, Saturations)
#     @tullio kr[ph, i] = relative_permeability(get_sat(kr_def, ph, i), Saturations[ph, i])
# end
