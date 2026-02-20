import{_ as s,c as n,a5 as e,o as p}from"./chunks/framework.DV8_rcIL.js";const l="/JutulDarcy.jl/previews/PR80/assets/pnednqb.DZsiuHru.jpeg",m=JSON.parse('{"title":"OLYMPUS_1 model","description":"","frontmatter":{},"headers":[],"relativePath":"examples/validation_olympus_1.md","filePath":"examples/validation_olympus_1.md","lastUpdated":null}'),t={name:"examples/validation_olympus_1.md"};function i(o,a,r,c,d,u){return p(),n("div",null,a[0]||(a[0]=[e(`<h1 id="OLYMPUS_1-model" tabindex="-1">OLYMPUS_1 model <a class="header-anchor" href="#OLYMPUS_1-model" aria-label="Permalink to &quot;OLYMPUS_1 model {#OLYMPUS_1-model}&quot;">​</a></h1><p>Model from the <a href="https://www.isapp2.com/optimization-challenge/reservoir-model-description.html" target="_blank" rel="noreferrer">ISAPP Optimization challenge</a></p><p>Two-phase complex corner-point model with primary and secondary production.</p><p>For more details, see [<a href="/JutulDarcy.jl/previews/PR80/extras/refs#olympus">10</a>]</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>using Jutul, JutulDarcy, GLMakie, DelimitedFiles, HYPRE</span></span>
<span class="line"><span>olympus_dir = JutulDarcy.GeoEnergyIO.test_input_file_path(&quot;OLYMPUS_1&quot;)</span></span>
<span class="line"><span>case = setup_case_from_data_file(joinpath(olympus_dir, &quot;OLYMPUS_1.DATA&quot;), backend = :csr)</span></span>
<span class="line"><span>ws, states = simulate_reservoir(case, output_substates = true)</span></span></code></pre></div><h2 id="Plot-the-reservoir" tabindex="-1">Plot the reservoir <a class="header-anchor" href="#Plot-the-reservoir" aria-label="Permalink to &quot;Plot the reservoir {#Plot-the-reservoir}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">model, key </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :porosity</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="`+l+`" alt=""></p><h2 id="Plot-the-saturations" tabindex="-1">Plot the saturations <a class="header-anchor" href="#Plot-the-saturations" aria-label="Permalink to &quot;Plot the saturations {#Plot-the-saturations}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>plot_reservoir(case.model, states, step = 10, key = :Saturations)</span></span></code></pre></div><h2 id="Load-reference-and-set-up-plotting" tabindex="-1">Load reference and set up plotting <a class="header-anchor" href="#Load-reference-and-set-up-plotting" aria-label="Permalink to &quot;Load reference and set up plotting {#Load-reference-and-set-up-plotting}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>csv_path = joinpath(olympus_dir, &quot;REFERENCE.CSV&quot;)</span></span>
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
<span class="line"><span>function plot_well_comparison(response, well_names, reponse_name = &quot;$response&quot;; ylims = missing)</span></span>
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
<span class="line"><span></span></span>
<span class="line"><span></span></span>
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
<span class="line"><span>        if !ismissing(ylims)</span></span>
<span class="line"><span>            ylims!(ax, ylims)</span></span>
<span class="line"><span>        end</span></span>
<span class="line"><span>    end</span></span>
<span class="line"><span>    l1 = LineElement(color = :black, linestyle = nothing)</span></span>
<span class="line"><span>    l2 = LineElement(color = :black, linestyle = :dash)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    Legend(fig[1:3, 2], linehandles, linelabels, nbanks = 3)</span></span>
<span class="line"><span>    Legend(fig[4, 2], [l1, l2], [&quot;JutulDarcy.jl&quot;, &quot;OPM Flow&quot;])</span></span>
<span class="line"><span>    fig</span></span>
<span class="line"><span>end</span></span></code></pre></div><h2 id="Plot-water-production-rates" tabindex="-1">Plot water production rates <a class="header-anchor" href="#Plot-water-production-rates" aria-label="Permalink to &quot;Plot water production rates {#Plot-water-production-rates}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>plot_well_comparison(:wrat, prod, &quot;Water surface production rate&quot;)</span></span></code></pre></div><h2 id="Plot-oil-production-rates" tabindex="-1">Plot oil production rates <a class="header-anchor" href="#Plot-oil-production-rates" aria-label="Permalink to &quot;Plot oil production rates {#Plot-oil-production-rates}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>plot_well_comparison(:orat, prod, &quot;Oil surface production rate&quot;, ylims = (0, 5))</span></span></code></pre></div><h2 id="Plot-water-injection-rates" tabindex="-1">Plot water injection rates <a class="header-anchor" href="#Plot-water-injection-rates" aria-label="Permalink to &quot;Plot water injection rates {#Plot-water-injection-rates}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>plot_well_comparison(:wrat, inj, &quot;Water injection rate&quot;)</span></span></code></pre></div><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/validation_olympus_1.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/validation_olympus_1.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,22)]))}const g=s(t,[["render",i]]);export{m as __pageData,g as default};
