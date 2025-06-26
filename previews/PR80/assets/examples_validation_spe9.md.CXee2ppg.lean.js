import{_ as s,c as n,a5 as e,o as p}from"./chunks/framework.DV8_rcIL.js";const h=JSON.parse('{"title":"The SPE9 model","description":"","frontmatter":{},"headers":[],"relativePath":"examples/validation_spe9.md","filePath":"examples/validation_spe9.md","lastUpdated":null}'),l={name:"examples/validation_spe9.md"};function t(i,a,r,o,c,u){return p(),n("div",null,a[0]||(a[0]=[e(`<h1 id="The-SPE9-model" tabindex="-1">The SPE9 model <a class="header-anchor" href="#The-SPE9-model" aria-label="Permalink to &quot;The SPE9 model {#The-SPE9-model}&quot;">​</a></h1><p><a href="http://dx.doi.org/10.2118/29110-MS" target="_blank" rel="noreferrer">Killough, J. E. 1995. Ninth SPE comparative solution project: A reexamination of black-oil simulation. In SPE Reservoir Simulation Symposium, 12-15 February 1995, San Antonio, Texas. SPE 29110-MS</a> [<a href="/JutulDarcy.jl/previews/PR80/extras/refs#spe9">8</a>]</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>using Jutul, JutulDarcy, GLMakie, DelimitedFiles</span></span>
<span class="line"><span>spe9_dir = JutulDarcy.GeoEnergyIO.test_input_file_path(&quot;SPE9&quot;)</span></span>
<span class="line"><span>case = setup_case_from_data_file(joinpath(spe9_dir, &quot;SPE9.DATA&quot;))</span></span>
<span class="line"><span>ws, states = simulate_reservoir(case, output_substates = true)</span></span>
<span class="line"><span>#</span></span>
<span class="line"><span>csv_path = joinpath(spe9_dir, &quot;REFERENCE.CSV&quot;)</span></span>
<span class="line"><span>data, header = readdlm(csv_path, &#39;,&#39;, header = true)</span></span>
<span class="line"><span>#</span></span>
<span class="line highlighted"><span>time_ref = data[:, 1]</span></span>
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
<span class="line"><span>    qts = well[:wrat] + well[:orat] + well[:grat]</span></span>
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
<span class="line"><span>    elseif response == :grat</span></span>
<span class="line"><span>        ys = si_unit(:day)/1e6</span></span>
<span class="line"><span>        yl = &quot;Surface gas rate / 10⁶ m³/day&quot;</span></span>
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
<span class="line"><span>        grat_ref = data[:, findfirst(x -&gt; x == &quot;$well_name:grat&quot;, vec(header))]</span></span>
<span class="line"><span>        orat_ref = data[:, findfirst(x -&gt; x == &quot;$well_name:orat&quot;, vec(header))]</span></span>
<span class="line"><span>        wrat_ref = data[:, findfirst(x -&gt; x == &quot;$well_name:wrat&quot;, vec(header))]</span></span>
<span class="line"><span>        tot_rate_ref = grat_ref + orat_ref + wrat_ref</span></span>
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
<span class="line"><span>end</span></span></code></pre></div><h2 id="Producer-BHP" tabindex="-1">Producer BHP <a class="header-anchor" href="#Producer-BHP" aria-label="Permalink to &quot;Producer BHP {#Producer-BHP}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot_well_comparison(:bhp, prod, &quot;Bottom hole pressure&quot;)</span></span></code></pre></div><h2 id="Producer-water-rate" tabindex="-1">Producer water rate <a class="header-anchor" href="#Producer-water-rate" aria-label="Permalink to &quot;Producer water rate {#Producer-water-rate}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot_well_comparison(:wrat, prod, &quot;Water surface rate&quot;)</span></span></code></pre></div><h2 id="Injector-water-rate" tabindex="-1">Injector water rate <a class="header-anchor" href="#Injector-water-rate" aria-label="Permalink to &quot;Injector water rate {#Injector-water-rate}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot_well_comparison(:wrat, inj, &quot;Water surface rate&quot;)</span></span></code></pre></div><h2 id="Oil-produciton-rate" tabindex="-1">Oil produciton rate <a class="header-anchor" href="#Oil-produciton-rate" aria-label="Permalink to &quot;Oil produciton rate {#Oil-produciton-rate}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot_well_comparison(:orat, prod, &quot;Oil surface rate&quot;)</span></span></code></pre></div><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/validation_spe9.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/validation_spe9.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,15)]))}const m=s(l,[["render",t]]);export{h as __pageData,m as default};
