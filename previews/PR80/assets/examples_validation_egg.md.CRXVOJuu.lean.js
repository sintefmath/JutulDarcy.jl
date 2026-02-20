import{_ as s,c as n,a5 as e,o as l}from"./chunks/framework.DV8_rcIL.js";const h=JSON.parse('{"title":"Validation: The Egg model (oil-water compressible)","description":"","frontmatter":{},"headers":[],"relativePath":"examples/validation_egg.md","filePath":"examples/validation_egg.md","lastUpdated":null}'),p={name:"examples/validation_egg.md"};function t(i,a,o,r,c,d){return l(),n("div",null,a[0]||(a[0]=[e(`<h1 id="Validation:-The-Egg-model-(oil-water-compressible)" tabindex="-1">Validation: The Egg model (oil-water compressible) <a class="header-anchor" href="#Validation:-The-Egg-model-(oil-water-compressible)" aria-label="Permalink to &quot;Validation: The Egg model (oil-water compressible) {#Validation:-The-Egg-model-(oil-water-compressible)}&quot;">​</a></h1><p>A two-phase model that is taken from the first member of the EGG ensemble. The model is a synthetic case with channelized permeability and water injection with fixed controls. For more details, see the paper where the ensemble is introduced:</p><p><a href="https://doi.org/10.1002/gdj3.21" target="_blank" rel="noreferrer">Jansen, Jan-Dirk, et al. &quot;The egg model–a geological ensemble for reservoir simulation.&quot; Geoscience Data Journal 1.2 (2014): 192-195.</a></p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>using Jutul, JutulDarcy, GLMakie, DelimitedFiles, HYPRE</span></span>
<span class="line"><span>egg_dir = JutulDarcy.GeoEnergyIO.test_input_file_path(&quot;EGG&quot;)</span></span>
<span class="line"><span>case = setup_case_from_data_file(joinpath(egg_dir, &quot;EGG.DATA&quot;))</span></span>
<span class="line"><span>ws, states = simulate_reservoir(case, output_substates = true)</span></span></code></pre></div><h2 id="Plot-the-reservoir-solution" tabindex="-1">Plot the reservoir solution <a class="header-anchor" href="#Plot-the-reservoir-solution" aria-label="Permalink to &quot;Plot the reservoir solution {#Plot-the-reservoir-solution}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot_reservoir(case.model, states, step = 135, key = :Saturations)</span></span></code></pre></div><h2 id="Load-reference-solution-(OPM-Flow)" tabindex="-1">Load reference solution (OPM Flow) <a class="header-anchor" href="#Load-reference-solution-(OPM-Flow)" aria-label="Permalink to &quot;Load reference solution (OPM Flow) {#Load-reference-solution-(OPM-Flow)}&quot;">​</a></h2><p>We load a CSV file with the reference solution and set up plotting</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>csv_path = joinpath(egg_dir, &quot;REFERENCE.CSV&quot;)</span></span>
<span class="line"><span>data, header = readdlm(csv_path, &#39;,&#39;, header = true)</span></span>
<span class="line"><span>time_ref = data[:, 1]</span></span>
<span class="line"><span>time_jutul = deepcopy(ws.time)</span></span>
<span class="line"><span>wells = deepcopy(ws.wells)</span></span>
<span class="line"><span>wnames = collect(keys(wells))</span></span>
<span class="line"><span>nw = length(wnames)</span></span>
<span class="line"><span>day = si_unit(:day)</span></span>
<span class="line"><span>cmap = :tableau_hue_circle</span></span>
<span class="line"><span></span></span>
<span class="line"><span>inj = Symbol[]</span></span>
<span class="line"><span>prod = Symbol[]</span></span>
<span class="line"><span>for (wellname, well) in pairs(wells)</span></span>
<span class="line"><span>    qts = well[:wrat] + well[:orat]</span></span>
<span class="line"><span>    if sum(qts) &gt; 0</span></span>
<span class="line"><span>        push!(inj, wellname)</span></span>
<span class="line"><span>    else</span></span>
<span class="line"><span>        push!(prod, wellname)</span></span>
<span class="line"><span>    end</span></span>
<span class="line"><span>end</span></span>
<span class="line"><span></span></span>
<span class="line"><span>function plot_well_comparison(response, well_names, reponse_name = &quot;$response&quot;)</span></span>
<span class="line"><span>    fig = Figure(size = (1000, 400))</span></span>
<span class="line"><span>    if response == :bhp</span></span>
<span class="line"><span>        ys = 1/si_unit(:bar)</span></span>
<span class="line"><span>        yl = &quot;Bottom hole pressure / Bar&quot;</span></span>
<span class="line"><span>    elseif response == :wrat</span></span>
<span class="line"><span>        ys = si_unit(:day)</span></span>
<span class="line"><span>        yl = &quot;Surface water rate / m³/day&quot;</span></span>
<span class="line"><span>    elseif response == :orat</span></span>
<span class="line"><span>        ys = si_unit(:day)/(1000*si_unit(:stb))</span></span>
<span class="line"><span>        yl = &quot;Surface oil rate / 10³ stb/day&quot;</span></span>
<span class="line"><span>    else</span></span>
<span class="line"><span>        error(&quot;$response not ready.&quot;)</span></span>
<span class="line"><span>    end</span></span>
<span class="line"><span>    welltypes = []</span></span>
<span class="line"><span>    ax = Axis(fig[1:4, 1], xlabel = &quot;Time / days&quot;, ylabel = yl)</span></span>
<span class="line"><span>    i = 1</span></span>
<span class="line"><span>    linehandles = []</span></span>
<span class="line"><span>    linelabels = []</span></span>
<span class="line"><span>    for well_name in well_names</span></span>
<span class="line"><span>        well = wells[well_name]</span></span>
<span class="line"><span>        label_in_csv = &quot;$well_name:$response&quot;</span></span>
<span class="line"><span>        ref_pos = findfirst(x -&gt; x == label_in_csv, vec(header))</span></span>
<span class="line"><span>        qoi = copy(well[response]).*ys</span></span>
<span class="line"><span>        qoi_ref = data[:, ref_pos].*ys</span></span>
<span class="line"><span>        tot_rate = copy(well[:rate])</span></span>
<span class="line"><span>        @. qoi[tot_rate == 0] = NaN</span></span>
<span class="line"><span>        orat_ref = data[:, findfirst(x -&gt; x == &quot;$well_name:orat&quot;, vec(header))]</span></span>
<span class="line"><span>        wrat_ref = data[:, findfirst(x -&gt; x == &quot;$well_name:wrat&quot;, vec(header))]</span></span>
<span class="line"><span>        tot_rate_ref = orat_ref + wrat_ref</span></span>
<span class="line"><span>        @. qoi_ref[tot_rate_ref == 0] = NaN</span></span>
<span class="line"><span>        crange = (1, max(length(well_names), 2))</span></span>
<span class="line"><span>        lh = lines!(ax, time_jutul./day, abs.(qoi),</span></span>
<span class="line"><span>            color = i,</span></span>
<span class="line"><span>            colorrange = crange,</span></span>
<span class="line"><span>            label = &quot;$well_name&quot;, colormap = cmap</span></span>
<span class="line"><span>        )</span></span>
<span class="line"><span>        push!(linehandles, lh)</span></span>
<span class="line"><span>        push!(linelabels, &quot;$well_name&quot;)</span></span>
<span class="line"><span>        lines!(ax, time_ref./day, abs.(qoi_ref),</span></span>
<span class="line"><span>            color = i,</span></span>
<span class="line"><span>            colorrange = crange,</span></span>
<span class="line"><span>            linestyle = :dash,</span></span>
<span class="line"><span>            colormap = cmap</span></span>
<span class="line"><span>        )</span></span>
<span class="line"><span>        i += 1</span></span>
<span class="line"><span>    end</span></span>
<span class="line"><span>    l1 = LineElement(color = :black, linestyle = nothing)</span></span>
<span class="line"><span>    l2 = LineElement(color = :black, linestyle = :dash)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    Legend(fig[1:3, 2], linehandles, linelabels, nbanks = 3)</span></span>
<span class="line"><span>    Legend(fig[4, 2], [l1, l2], [&quot;JutulDarcy.jl&quot;, &quot;E100&quot;])</span></span>
<span class="line"><span>    fig</span></span>
<span class="line"><span>end</span></span></code></pre></div><h2 id="Well-responses-and-comparison" tabindex="-1">Well responses and comparison <a class="header-anchor" href="#Well-responses-and-comparison" aria-label="Permalink to &quot;Well responses and comparison {#Well-responses-and-comparison}&quot;">​</a></h2><p>As the case is a two-phase model with water injection, we limit the results to plots of the producer water and oil rates.</p><h3 id="Water-production-rates" tabindex="-1">Water production rates <a class="header-anchor" href="#Water-production-rates" aria-label="Permalink to &quot;Water production rates {#Water-production-rates}&quot;">​</a></h3><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot_well_comparison(:wrat, prod, &quot;Producer water surface rate&quot;)</span></span></code></pre></div><h3 id="Oil-production-rates" tabindex="-1">Oil production rates <a class="header-anchor" href="#Oil-production-rates" aria-label="Permalink to &quot;Oil production rates {#Oil-production-rates}&quot;">​</a></h3><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot_well_comparison(:orat, prod, &quot;Producer oil surface rate&quot;)</span></span></code></pre></div><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/validation_egg.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/validation_egg.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,19)]))}const m=s(p,[["render",t]]);export{h as __pageData,m as default};
