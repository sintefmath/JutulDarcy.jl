using Jutul, JutulDarcy, TOML
include(joinpath(@__DIR__, "refine_layers.jl"))
function same_entry(a,b)
    if a isa AbstractArray && b isa AbstractArray
        return size(a) == size(b) && all(same_entry.(a,b))
    elseif a isa AbstractDict && b isa AbstractDict
        return keys(a) == keys(b) && all(same_entry(a[k],b[k]) for k in keys(a))
    end
    return isequal(a,b) || (a isa Number && b isa Number && isapprox(a,b; rtol=5e-14, atol=0))
end
results = Any[]
for (name, factor, stem) in (("olympus",4,"OLYMPUS_1"),("norne",8,"NORNE_NOHYST"))
    path = JutulDarcy.GeoEnergyIO.test_input_file_path(stem, "$stem.DATA")
    expected = refine_layers(JutulDarcy.GeoEnergyIO.parse_data_file(path), factor)
    exported = joinpath(@__DIR__, "results", "opm", "$(name)_z$(factor)", uppercase("$(name)_Z$(factor)")*".DATA")
    actual = JutulDarcy.GeoEnergyIO.parse_data_file(exported)
    checked = String[]
    for section in ("RUNSPEC","GRID","PROPS","REGIONS","SOLUTION")
        for (key, a) in expected[section]
            a isa AbstractArray || continue
            key in ("COORDSYS", "WELLDIMS", "EQLDIMS") && continue
            haskey(actual[section], key) || error("Missing $section/$key")
            b = actual[section][key]
            @assert size(a) == size(b) (section,key,size(a),size(b))
            @assert same_entry(a,b) (section,key)
            push!(checked,"$section/$key")
        end
    end
    for key in ("FAULTS", "MINPV")
        if haskey(expected["GRID"],key)
            @assert same_entry(expected["GRID"][key],actual["GRID"][key]) key
        end
    end
    target_steps = actual["SCHEDULE"]["STEPS"]
    source_steps = expected["SCHEDULE"]["STEPS"]
    completions = 0
    for (i, b) in enumerate(target_steps)
        isempty(b) && continue
        a = source_steps[i]
        for key in ("COMPDAT","WELOPEN","DATES","TSTEP","WELSPECS","WCONPROD","WCONINJE","WCONHIST","WCONINJH","WELTARG","WEFAC","GRUPTREE","WGRUPCON","GCONPROD","GCONINJE")
            av, bv = get(a,key,[]), get(b,key,[])
            @assert length(av) == length(bv) (i,key,length(av),length(bv))
            for (left,right) in zip(av,bv)
                if left isa AbstractArray
                    @assert length(left) == length(right) (i,key)
                    @assert all(same_entry.(left,right)) (i,key,left,right)
                else
                    @assert same_entry(left,right) (i,key,left,right)
                end
            end
        end
        completions += length(get(b,"COMPDAT",[]))
    end
    row = Dict("case"=>name,"factor"=>factor,"fields"=>checked,"schedule_blocks"=>length(target_steps),"completions"=>completions)
    push!(results,row); println(row); flush(stdout)
end
open(joinpath(@__DIR__,"results","opm_deck_checks.toml"),"w") do io
    TOML.print(io,Dict("results"=>results))
end
