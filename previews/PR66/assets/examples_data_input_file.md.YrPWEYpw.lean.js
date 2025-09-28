import{_ as a,c as n,a5 as t,o as e}from"./chunks/framework.zT1yXYO7.js";const i="/JutulDarcy.jl/previews/PR66/assets/eoockhv.B-YoMvy0.jpeg",p="/JutulDarcy.jl/previews/PR66/assets/mzzvzjy.BE1NmEEx.jpeg",k=JSON.parse('{"title":"Simulating Eclipse/DATA input files","description":"","frontmatter":{},"headers":[],"relativePath":"examples/data_input_file.md","filePath":"examples/data_input_file.md","lastUpdated":null}'),l={name:"examples/data_input_file.md"};function o(r,s,c,h,u,d){return e(),n("div",null,s[0]||(s[0]=[t(`<h1 id="Simulating-Eclipse/DATA-input-files" tabindex="-1">Simulating Eclipse/DATA input files <a class="header-anchor" href="#Simulating-Eclipse/DATA-input-files" aria-label="Permalink to &quot;Simulating Eclipse/DATA input files {#Simulating-Eclipse/DATA-input-files}&quot;">​</a></h1><p>The DATA format is commonly used in reservoir simulation. JutulDarcy can set up cases on this format and includes a fully featured grid builder for corner-point grids. Once a case has been set up, it uses the same types as a regular JutulDarcy simulation, allowing modification and use of the case in differentiable workflows.</p><p>We begin by loading the SPE9 dataset via the GeoEnergyIO package.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">pth </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">GeoEnergyIO</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">test_input_file_path</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;SPE9&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;SPE9.DATA&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>&quot;/home/runner/.julia/artifacts/c1063fdc96b21cbc18ea8a3bd1a1281aa04ffa3f/SPE9.DATA&quot;</span></span></code></pre></div><h2 id="Set-up-and-run-a-simulation" tabindex="-1">Set up and run a simulation <a class="header-anchor" href="#Set-up-and-run-a-simulation" aria-label="Permalink to &quot;Set up and run a simulation {#Set-up-and-run-a-simulation}&quot;">​</a></h2><p>If we do not need the case, we could also have done: ws, states = simulate_data_file(pth)</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">case </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_case_from_data_file</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(pth)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ws, states </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>ReservoirSimResult with 90 entries:</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  wells (26 present):</span></span>
<span class="line"><span>    :INJE1</span></span>
<span class="line"><span>    :PRODU25</span></span>
<span class="line"><span>    :PRODU16</span></span>
<span class="line"><span>    :PRODU20</span></span>
<span class="line"><span>    :PRODU10</span></span>
<span class="line"><span>    :PRODU12</span></span>
<span class="line"><span>    :PRODU22</span></span>
<span class="line"><span>    :PRODU14</span></span>
<span class="line"><span>    :PRODU6</span></span>
<span class="line"><span>    :PRODU7</span></span>
<span class="line"><span>    :PRODU9</span></span>
<span class="line"><span>    :PRODU24</span></span>
<span class="line"><span>    :PRODU3</span></span>
<span class="line"><span>    :PRODU23</span></span>
<span class="line"><span>    :PRODU5</span></span>
<span class="line"><span>    :PRODU11</span></span>
<span class="line"><span>    :PRODU17</span></span>
<span class="line"><span>    :PRODU4</span></span>
<span class="line"><span>    :PRODU15</span></span>
<span class="line"><span>    :PRODU2</span></span>
<span class="line"><span>    :PRODU19</span></span>
<span class="line"><span>    :PRODU21</span></span>
<span class="line"><span>    :PRODU13</span></span>
<span class="line"><span>    :PRODU8</span></span>
<span class="line"><span>    :PRODU26</span></span>
<span class="line"><span>    :PRODU18</span></span>
<span class="line"><span>    Results per well:</span></span>
<span class="line"><span>       :wrat =&gt; Vector{Float64} of size (90,)</span></span>
<span class="line"><span>       :Aqueous_mass_rate =&gt; Vector{Float64} of size (90,)</span></span>
<span class="line"><span>       :orat =&gt; Vector{Float64} of size (90,)</span></span>
<span class="line"><span>       :bhp =&gt; Vector{Float64} of size (90,)</span></span>
<span class="line"><span>       :lrat =&gt; Vector{Float64} of size (90,)</span></span>
<span class="line"><span>       :mass_rate =&gt; Vector{Float64} of size (90,)</span></span>
<span class="line"><span>       :rate =&gt; Vector{Float64} of size (90,)</span></span>
<span class="line"><span>       :Vapor_mass_rate =&gt; Vector{Float64} of size (90,)</span></span>
<span class="line"><span>       :control =&gt; Vector{Symbol} of size (90,)</span></span>
<span class="line"><span>       :Liquid_mass_rate =&gt; Vector{Float64} of size (90,)</span></span>
<span class="line"><span>       :grat =&gt; Vector{Float64} of size (90,)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  states (Vector with 90 entries, reservoir variables for each state)</span></span>
<span class="line"><span>    :BlackOilUnknown =&gt; Vector{BlackOilX{Float64}} of size (9000,)</span></span>
<span class="line"><span>    :Saturations =&gt; Matrix{Float64} of size (3, 9000)</span></span>
<span class="line"><span>    :Pressure =&gt; Vector{Float64} of size (9000,)</span></span>
<span class="line"><span>    :Rs =&gt; Vector{Float64} of size (9000,)</span></span>
<span class="line"><span>    :ImmiscibleSaturation =&gt; Vector{Float64} of size (9000,)</span></span>
<span class="line"><span>    :TotalMasses =&gt; Matrix{Float64} of size (3, 9000)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  time (report time for each state)</span></span>
<span class="line"><span>     Vector{Float64} of length 90</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  result (extended states, reports)</span></span>
<span class="line"><span>     SimResult with 90 entries</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  extra</span></span>
<span class="line"><span>     Dict{Any, Any} with keys :simulator, :config</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Completed at Oct. 01 2024 13:10 after 12 seconds, 552 milliseconds, 585 microseconds.</span></span></code></pre></div><h2 id="Show-the-input-data" tabindex="-1">Show the input data <a class="header-anchor" href="#Show-the-input-data" aria-label="Permalink to &quot;Show the input data {#Show-the-input-data}&quot;">​</a></h2><p>The input data takes the form of a Dict:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">case</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">input_data</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Dict{String, Any} with 6 entries:</span></span>
<span class="line"><span>  &quot;RUNSPEC&quot;  =&gt; OrderedDict{String, Any}(&quot;TITLE&quot;=&gt;&quot;SPE 9&quot;, &quot;DIMENS&quot;=&gt;[24, 25, 1…</span></span>
<span class="line"><span>  &quot;GRID&quot;     =&gt; OrderedDict{String, Any}(&quot;cartDims&quot;=&gt;(24, 25, 15), &quot;CURRENT_BOX…</span></span>
<span class="line"><span>  &quot;PROPS&quot;    =&gt; OrderedDict{String, Any}(&quot;PVTW&quot;=&gt;Any[[2.48211e7, 1.0034, 4.3511…</span></span>
<span class="line"><span>  &quot;SUMMARY&quot;  =&gt; OrderedDict{String, Any}()</span></span>
<span class="line"><span>  &quot;SCHEDULE&quot; =&gt; Dict{String, Any}(&quot;STEPS&quot;=&gt;OrderedDict{String, Any}[OrderedDict…</span></span>
<span class="line"><span>  &quot;SOLUTION&quot; =&gt; OrderedDict{String, Any}(&quot;EQUIL&quot;=&gt;Any[[2753.87, 2.48211e7, 3032…</span></span></code></pre></div><p>We can also examine the RUNSPEC section</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">case</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">input_data[</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;RUNSPEC&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>OrderedCollections.OrderedDict{String, Any} with 13 entries:</span></span>
<span class="line"><span>  &quot;TITLE&quot;    =&gt; &quot;SPE 9&quot;</span></span>
<span class="line"><span>  &quot;DIMENS&quot;   =&gt; [24, 25, 15]</span></span>
<span class="line"><span>  &quot;OIL&quot;      =&gt; true</span></span>
<span class="line"><span>  &quot;WATER&quot;    =&gt; true</span></span>
<span class="line"><span>  &quot;GAS&quot;      =&gt; true</span></span>
<span class="line"><span>  &quot;DISGAS&quot;   =&gt; true</span></span>
<span class="line"><span>  &quot;FIELD&quot;    =&gt; true</span></span>
<span class="line"><span>  &quot;START&quot;    =&gt; DateTime(&quot;2015-01-01T00:00:00&quot;)</span></span>
<span class="line"><span>  &quot;WELLDIMS&quot; =&gt; [26, 5, 1, 26, 5, 10, 5, 4, 3, 0, 1, 1, 10, 201]</span></span>
<span class="line"><span>  &quot;TABDIMS&quot;  =&gt; [1, 1, 40, 20, 1, 20, 20, 1, 1, -1  …  -1, 10, 10, 10, -1, 5, 5…</span></span>
<span class="line"><span>  &quot;EQLDIMS&quot;  =&gt; [1, 100, 50, 1, 50]</span></span>
<span class="line"><span>  &quot;UNIFIN&quot;   =&gt; true</span></span>
<span class="line"><span>  &quot;UNIFOUT&quot;  =&gt; true</span></span></code></pre></div><h2 id="Plot-the-simulation-model" tabindex="-1">Plot the simulation model <a class="header-anchor" href="#Plot-the-simulation-model" aria-label="Permalink to &quot;Plot the simulation model {#Plot-the-simulation-model}&quot;">​</a></h2><p>These plot are interactive when run outside of the documentations.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GLMakie</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">model, states)</span></span></code></pre></div><p><img src="`+i+'" alt=""></p><h2 id="Plot-the-well-responses" tabindex="-1">Plot the well responses <a class="header-anchor" href="#Plot-the-well-responses" aria-label="Permalink to &quot;Plot the well responses {#Plot-the-well-responses}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_well_results</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ws)</span></span></code></pre></div><p><img src="'+p+'" alt=""></p><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/data_input_file.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/data_input_file.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',27)]))}const E=a(l,[["render",o]]);export{k as __pageData,E as default};
