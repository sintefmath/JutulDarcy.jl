import{_ as l,c as e,j as a,a as n,a5 as t,o as p}from"./chunks/framework.DV8_rcIL.js";const v=JSON.parse('{"title":"Validation of equation-of-state compositional simulator","description":"","frontmatter":{},"headers":[],"relativePath":"examples/validation_compositional.md","filePath":"examples/validation_compositional.md","lastUpdated":null}'),i={name:"examples/validation_compositional.md"},o={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},r={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"0.988ex",height:"1.065ex",role:"img",focusable:"false",viewBox:"0 -320.9 436.6 470.9","aria-hidden":"true"},c={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},u={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"0.988ex",height:"1.065ex",role:"img",focusable:"false",viewBox:"0 -320.9 436.6 470.9","aria-hidden":"true"};function m(d,s,h,g,x,f){return p(),e("div",null,[s[9]||(s[9]=a("h1",{id:"Validation-of-equation-of-state-compositional-simulator",tabindex:"-1"},[n("Validation of equation-of-state compositional simulator "),a("a",{class:"header-anchor",href:"#Validation-of-equation-of-state-compositional-simulator","aria-label":'Permalink to "Validation of equation-of-state compositional simulator {#Validation-of-equation-of-state-compositional-simulator}"'},"​")],-1)),s[10]||(s[10]=a("p",null,"This example solves a 1D two-phase, three component miscible displacement problem and compares against existing simulators (E300, AD-GPRS) to verify correctness.",-1)),s[11]||(s[11]=a("p",null,"The case is loaded from an input file that can be run in other simulators. For convenience, we provide solutions from the other simulators as a binary file to perform a comparison without having to run and convert results from other the simulators.",-1)),a("p",null,[s[4]||(s[4]=n("This case is a small compositional problem inspired by the examples in Voskov et al (JPSE, 2012). A 1D reservoir of 1,000 meters length is discretized into 1,000 cells. The model initially contains a mixture made up of 0.6 parts C10, 0.1 parts CO2, and 0.3 parts C1 by moles at 150 degrees C and 75 bar pressure. Wells are placed in the leftmost and rightmost cells of the domain, with the leftmost well injecting pure CO")),a("mjx-container",o,[(p(),e("svg",r,s[0]||(s[0]=[t('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="msub"><g data-mml-node="mi"></g><g data-mml-node="mn" transform="translate(33,-150) scale(0.707)"><path data-c="32" d="M109 429Q82 429 66 447T50 491Q50 562 103 614T235 666Q326 666 387 610T449 465Q449 422 429 383T381 315T301 241Q265 210 201 149L142 93L218 92Q375 92 385 97Q392 99 409 186V189H449V186Q448 183 436 95T421 3V0H50V19V31Q50 38 56 46T86 81Q115 113 136 137Q145 147 170 174T204 211T233 244T261 278T284 308T305 340T320 369T333 401T340 431T343 464Q343 527 309 573T212 619Q179 619 154 602T119 569T109 550Q109 549 114 549Q132 549 151 535T170 489Q170 464 154 447T109 429Z" style="stroke-width:3;"></path></g></g></g></g>',1)]))),s[1]||(s[1]=a("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[a("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[a("msub",null,[a("mi"),a("mn",null,"2")])])],-1))]),s[5]||(s[5]=n(" at a fixed bottom-hole pressure of 100 bar and the other well producing at 50 bar. The model is isothermal and contains a phase transition from the initial two-phase mixture to single-phase gas as injected CO")),a("mjx-container",c,[(p(),e("svg",u,s[2]||(s[2]=[t('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="msub"><g data-mml-node="mi"></g><g data-mml-node="mn" transform="translate(33,-150) scale(0.707)"><path data-c="32" d="M109 429Q82 429 66 447T50 491Q50 562 103 614T235 666Q326 666 387 610T449 465Q449 422 429 383T381 315T301 241Q265 210 201 149L142 93L218 92Q375 92 385 97Q392 99 409 186V189H449V186Q448 183 436 95T421 3V0H50V19V31Q50 38 56 46T86 81Q115 113 136 137Q145 147 170 174T204 211T233 244T261 278T284 308T305 340T320 369T333 401T340 431T343 464Q343 527 309 573T212 619Q179 619 154 602T119 569T109 550Q109 549 114 549Q132 549 151 535T170 489Q170 464 154 447T109 429Z" style="stroke-width:3;"></path></g></g></g></g>',1)]))),s[3]||(s[3]=a("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[a("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[a("msub",null,[a("mi"),a("mn",null,"2")])])],-1))]),s[6]||(s[6]=n(" eventually displaces the resident fluids. For further details on this setup, see Møyner and Tchelepi (SPE J. 2018) [")),s[7]||(s[7]=a("a",{href:"/JutulDarcy.jl/previews/PR80/extras/refs#moyner_tchelepi_2018"},"9",-1)),s[8]||(s[8]=n("]."))]),s[12]||(s[12]=t(`<div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>using JutulDarcy</span></span>
<span class="line"><span>using Jutul</span></span>
<span class="line"><span>using GLMakie</span></span>
<span class="line"><span>dpth = JutulDarcy.GeoEnergyIO.test_input_file_path(&quot;SIMPLE_COMP&quot;)</span></span>
<span class="line"><span>data_path = joinpath(dpth, &quot;SIMPLE_COMP.DATA&quot;)</span></span>
<span class="line"><span>case = setup_case_from_data_file(data_path)</span></span>
<span class="line"><span>result = simulate_reservoir(case, info_level = -1)</span></span>
<span class="line"><span>ws, states = result;</span></span>
<span class="line"><span>nothing #hide</span></span></code></pre></div><h2 id="Plot-solutions-and-compare" tabindex="-1">Plot solutions and compare <a class="header-anchor" href="#Plot-solutions-and-compare" aria-label="Permalink to &quot;Plot solutions and compare {#Plot-solutions-and-compare}&quot;">​</a></h2><p>The 1D displacement can be plotted as a line plot. We pick a step midway through the simulation and plot compositions, saturations and pressure.</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>cmap = :tableau_hue_circle</span></span>
<span class="line"><span>ref_path = joinpath(dpth, &quot;reference.jld2&quot;)</span></span>
<span class="line"><span>ref = Jutul.JLD2.load(ref_path)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>step_to_plot = 250</span></span>
<span class="line"><span>fig = with_theme(theme_latexfonts()) do</span></span>
<span class="line"><span>    x = reservoir_domain(case)[:cell_centroids][1, :]</span></span>
<span class="line"><span>    mz = 3</span></span>
<span class="line"><span>    ix = step_to_plot</span></span>
<span class="line"><span>    mt = :circle</span></span>
<span class="line"><span>    fig = Figure(size = (800, 400))</span></span>
<span class="line"><span>    ax = Axis(fig[2, 1], xlabel = &quot;Cell center / m&quot;)</span></span>
<span class="line"><span>    cnames = [&quot;DECANE&quot;, &quot;CO2&quot;, &quot;METHANE&quot;]</span></span>
<span class="line"><span>    cnames = [&quot;C₁₀&quot;, &quot;CO₂&quot;, &quot;C₁&quot;]</span></span>
<span class="line"><span>    cnames = [L&quot;\\text{C}_{10}&quot;, L&quot;\\text{CO}_2&quot;, L&quot;\\text{C}_1&quot;]</span></span>
<span class="line"><span>    lineh = []</span></span>
<span class="line"><span>    lnames = []</span></span>
<span class="line"><span>    crange = (1, 4)</span></span>
<span class="line"><span>    for i in range(crange...)</span></span>
<span class="line"><span>        if i == 4</span></span>
<span class="line"><span>            cname = L&quot;\\text{S}_g&quot;</span></span>
<span class="line"><span>            gprs = missing</span></span>
<span class="line"><span>            ecl = ref[&quot;e300&quot;][ix][:Saturations][2, :]</span></span>
<span class="line"><span>            ju = states[ix][:Saturations][2, :]</span></span>
<span class="line"><span>        else</span></span>
<span class="line"><span>            @assert i &lt;= 4</span></span>
<span class="line"><span>            ecl = ref[&quot;e300&quot;][ix][:OverallMoleFractions][i, :]</span></span>
<span class="line"><span>            gprs = ref[&quot;adgprs&quot;][ix][:OverallMoleFractions][i, :]</span></span>
<span class="line"><span>            ju = states[ix][:OverallMoleFractions][i, :]</span></span>
<span class="line"><span>            cname = cnames[i]</span></span>
<span class="line"><span>        end</span></span>
<span class="line"><span>        h = lines!(ax, x, ju, colormap = cmap, color = i, colorrange = crange, label = cname)</span></span>
<span class="line"><span>        push!(lnames, cname)</span></span>
<span class="line"><span>        push!(lineh, h)</span></span>
<span class="line"><span>        scatter!(ax, x, ecl, markersize = mz, colormap = cmap, color = i, colorrange = crange)</span></span>
<span class="line"><span>        if !ismissing(gprs)</span></span>
<span class="line"><span>            lines!(ax, x, gprs, colormap = cmap, color = i, colorrange = crange, linestyle = :dash)</span></span>
<span class="line"><span>        end</span></span>
<span class="line"><span>    end</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    l_ju = LineElement(color = :black, linestyle = nothing)</span></span>
<span class="line"><span>    l_ecl = MarkerElement(color = :black, markersize = mz, marker = mt)</span></span>
<span class="line"><span>    l_gprs = LineElement(color = :black, linestyle = :dash)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    Legend(</span></span>
<span class="line"><span>        fig[1, 1],</span></span>
<span class="line"><span>        [[l_ju, l_ecl, l_gprs], lineh],</span></span>
<span class="line"><span>        [[L&quot;\\text{JutulDarcy}&quot;, L&quot;\\text{E300}&quot;, L&quot;\\text{AD-GPRS}&quot;], lnames],</span></span>
<span class="line"><span>        [&quot;Simulator&quot;, &quot;Result&quot;],</span></span>
<span class="line"><span>        tellwidth = false,</span></span>
<span class="line"><span>        orientation = :horizontal,</span></span>
<span class="line"><span>    )</span></span>
<span class="line"><span>    fig</span></span>
<span class="line"><span>end</span></span></code></pre></div><h2 id="Calculate-sensitivities" tabindex="-1">Calculate sensitivities <a class="header-anchor" href="#Calculate-sensitivities" aria-label="Permalink to &quot;Calculate sensitivities {#Calculate-sensitivities}&quot;">​</a></h2><p>We demonstrate how the parameter sensitivities of an objective function can be calculated for a compositional model.</p><p>The objective function is taken to be the average gas saturation at a specific report step that was plotted in the previous section.</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>import Statistics: mean</span></span>
<span class="line"><span>import JutulDarcy: reservoir_sensitivities</span></span>
<span class="line"><span>function objective_function(model, state, Δt, step_i, forces)</span></span>
<span class="line"><span>    if step_i != step_to_plot</span></span>
<span class="line"><span>        return 0.0</span></span>
<span class="line"><span>    end</span></span>
<span class="line"><span>    sg = @view state.Reservoir.Saturations[2, :]</span></span>
<span class="line"><span>    return mean(sg)</span></span>
<span class="line"><span>end</span></span>
<span class="line"><span>data_domain_with_gradients = reservoir_sensitivities(case, result, objective_function, include_parameters = true)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>fig = with_theme(theme_latexfonts()) do</span></span>
<span class="line"><span>    x = reservoir_domain(case)[:cell_centroids][1, :]</span></span>
<span class="line"><span>    mz = 3</span></span>
<span class="line"><span>    ix = step_to_plot</span></span>
<span class="line"><span>    cmap = :Dark2_5</span></span>
<span class="line"><span>    cmap = :Accent_4</span></span>
<span class="line"><span>    cmap = :Spectral_4</span></span>
<span class="line"><span>    cmap = :tableau_hue_circle</span></span>
<span class="line"><span>    mt = :circle</span></span>
<span class="line"><span>    fig = Figure(size = (800, 400))</span></span>
<span class="line"><span>    normalize = x -&gt; x./(maximum(x) - minimum(x))</span></span>
<span class="line"><span>    logscale = x -&gt; sign.(x).*log10.(abs.(x))</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    ∂T = data_domain_with_gradients[:temperature]</span></span>
<span class="line"><span>    ∂ϕ = data_domain_with_gradients[:porosity]</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    ax1 = Axis(fig[2, 1], yticklabelcolor = :blue, xlabel = &quot;Cell center / m&quot;)</span></span>
<span class="line"><span>    ax2 = Axis(fig[2, 1], yticklabelcolor = :red, yaxisposition = :right)</span></span>
<span class="line"><span>    hidespines!(ax2)</span></span>
<span class="line"><span>    hidexdecorations!(ax2)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    l1 = lines!(ax1, x, ∂T, label = &quot;Temperature&quot;, color = :blue)</span></span>
<span class="line"><span>    l2 = lines!(ax2, x, ∂ϕ, label = &quot;Porosity&quot;, color = :red)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    Legend(fig[1, 1], [l1, l2], [&quot;Temperature&quot;, &quot;Porosity&quot;], &quot;Parameter&quot;, tellwidth = false, orientation = :horizontal)</span></span>
<span class="line"><span>    fig</span></span>
<span class="line"><span>end</span></span>
<span class="line"><span>fig</span></span></code></pre></div><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/validation_compositional.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/validation_compositional.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,12))])}const _=l(i,[["render",m]]);export{v as __pageData,_ as default};
