# # Comparison between JutulDarcy.jl and MRST
# This example contains validation of JutulDarcy.jl against
# [MRST](https://www.mrst.no). In general, minor differences are observed. These
# can be traced back to a combination of different internal timestepping done by
# the simulators and that JutulDarcy by default uses a multisegment well
# formulation while MRST uses a standard instantaneous equilibrium model without
# well bore storage terms. These differences are most evident when simulators
# start up.
#
# These cases have been exported using the MRST `jutul` module which can export
# MRST or Eclipse-type of cases to a JutulDarcy-compatible input format. They
# can then be simulated using [`simulate_mrst_case`](@ref).
using JutulDarcy, Jutul
using GLMakie

# ## Define a few utilities for plotting the MRST results
# We are going to compare well responses against pre-computed results stored
# inside the JutulDarcy module.
function mrst_case_path(name)
    base_path, = splitdir(pathof(JutulDarcy))
    joinpath(base_path, "..", "test", "mrst", "$(name).mat")
end

function mrst_solution(result)
    return result.extra[:mrst]["extra"][1]["mrst_solution"]
end

function mrst_well_index(mrst_result, k)
    return findfirst(isequal("$k"), vec(mrst_result["names"]))
end

function get_mrst_comparison(wdata, ref, wname, t = :bhp)
    yscale = "m³/s"
    if t == :bhp
        tname = :bhp
        mname = "bhp"
        yscale = "Pa"
    elseif t == :qos
        tname = :orat
        mname = "qOs"
    elseif t == :qws
        tname = :wrat
        mname = "qWs"
    elseif t == :qgs
        tname = :grat
        mname = "qGs"
    else
        error("Not supported: $t")
    end
    jutul = wdata[tname]
    mrst = ref[mname][:, mrst_well_index(ref, wname)]

    return (jutul, mrst, tname, yscale)
end

function plot_comparison(wsol, ref, rep_t, t, wells_keys = keys(wsol.wells))
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "time (days)")
    l = ""
    yscale = ""
    T = rep_t./(3600*24.0)
    for (w, d) in wsol.wells
        if !(w in wells_keys)
            continue
        end
        jutul, mrst, l, yscale = get_mrst_comparison(d, ref, w, t)
        lines!(ax, T, abs.(jutul), label = "$w")
        scatter!(ax, T, abs.(mrst), markersize = 8)
    end
    axislegend()
    ax.ylabel[] = "$l ($yscale)"
    fig
end
#-
# ## The Egg model (oil-water compressible)
# A two-phase model that is taken from the first member of the EGG ensemble. For
# more details, see the paper where the ensemble is introduced:
#
# [Jansen, Jan-Dirk, et al. "The egg model–a geological ensemble for reservoir
# simulation." Geoscience Data Journal 1.2 (2014): 192-195.](https://doi.org/10.1002/gdj3.21)

# ### Simulate model
egg = simulate_mrst_case(mrst_case_path("egg"), info_level = -1)
wells = egg.wells
rep_t = egg.time
ref = mrst_solution(egg);
# ### Compare well responses
injectors = [:INJECT1, :INJECT2, :INJECT3, :INJECT4, :INJECT5, :INJECT6, :INJECT7]
producers = [:PROD1, :PROD2, :PROD3, :PROD4]
#-
# #### Bottom hole pressures
plot_comparison(wells, ref, rep_t, :bhp, injectors)
#-
# #### Oil rates
plot_comparison(wells, ref, rep_t, :qos, producers)
#-
# #### Water rates
plot_comparison(wells, ref, rep_t, :qws, producers)

# ## SPE1 (black oil, disgas)
# A shortened version of the SPE1 benchmark case.
#
# [Odeh, A.S. 1981. Comparison of Solutions to a Three-Dimensional Black-Oil
# Reservoir Simulation Problem. J Pet Technol 33 (1): 13–25.
# SPE-9723-PA](http://dx.doi.org/10.2118/9723-PA)
#
# For comparison against other simulators, see the equivialent JutulSPE1 example
# in the Jutul module for [MRST](https://www.mrst.no)
# ### Simulate model
spe1 = simulate_mrst_case(mrst_case_path("spe1"), info_level = -1)
wells = spe1.wells
rep_t = spe1.time
ref = mrst_solution(spe1);
# ### Compare well responses
# #### Bottom hole pressures
#-
plot_comparison(wells, ref, rep_t, :bhp)
# #### Gas rates
#-
plot_comparison(wells, ref, rep_t, :qgs, [:PRODUCER])

# ## SPE3 (black oil, vapoil)
# A black-oil variant of the SPE3 benchmark case.
#
# [Kenyon, D. "Third SPE comparative solution project: gas cycling of retrograde
# condensate reservoirs." Journal of Petroleum Technology 39.08 (1987):
# 981-997](http://dx.doi.org/10.2118/12278-PA)

# ### Simulate model
spe3 = simulate_mrst_case(mrst_case_path("spe3"), info_level = -1)
wells = spe3.wells
rep_t = spe3.time
ref = mrst_solution(spe3);
# ### Compare well responses
# #### Bottom hole pressures
plot_comparison(wells, ref, rep_t, :bhp, [:PRODUCER])
#-
#-
# #### Gas rates
plot_comparison(wells, ref, rep_t, :qgs, [:PRODUCER])
#-
#-
# #### Oil rates
plot_comparison(wells, ref, rep_t, :qos, [:PRODUCER])
#-
#-

# ## SPE9 (black oil, disgas)
# Example of the SPE9 model exported from MRST running in JutulDarcy.
#
#   [Killough, J. E. 1995. Ninth SPE comparative solution project: A
#   reexamination of black-oil simulation. In SPE Reservoir Simulation
#   Symposium,  12-15 February 1995, San Antonio, Texas. SPE 29110-MS]
#   (http://dx.doi.org/10.2118/29110-MS)
#
# For comparison against other simulators, see the equivialent JutulSPE9 example 
# in the Jutul module for [MRST](https://www.mrst.no)

# ### Simulate model
spe9 = simulate_mrst_case(mrst_case_path("spe9"), info_level = -1)
wells = spe9.wells
rep_t = spe9.time
ref = mrst_solution(spe9);
# ### Compare well responses
injectors = [:INJE1]
producers = [Symbol("PROD$i") for i in 1:25]

# #### Injector water rate
plot_comparison(wells, ref, rep_t, :qws, injectors)
#-
# #### Oil rates
plot_comparison(wells, ref, rep_t, :qos, producers)
#-
# #### Water rates
plot_comparison(wells, ref, rep_t, :qws, producers)
#-
# #### Bottom hole pressures
plot_comparison(wells, ref, rep_t, :bhp, producers)
#-
