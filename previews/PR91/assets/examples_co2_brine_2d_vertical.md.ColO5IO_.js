import{_ as a,c as n,ai as i,o as p}from"./chunks/framework.nF1RJhOO.js";const l="/JutulDarcy.jl/previews/PR91/assets/okymxsg.BwgSd9yE.jpeg",e="/JutulDarcy.jl/previews/PR91/assets/rgoncni.vwPmYP-E.jpeg",t="/JutulDarcy.jl/previews/PR91/assets/lenybjx.ArfBpk4g.jpeg",h="/JutulDarcy.jl/previews/PR91/assets/ticpxit.b3J6fV39.jpeg",b=JSON.parse('{"title":"Intro to compositional flow","description":"","frontmatter":{},"headers":[],"relativePath":"examples/co2_brine_2d_vertical.md","filePath":"examples/co2_brine_2d_vertical.md","lastUpdated":null}'),k={name:"examples/co2_brine_2d_vertical.md"};function r(c,s,o,E,d,g){return p(),n("div",null,s[0]||(s[0]=[i(`<h1 id="Intro-to-compositional-flow" tabindex="-1">Intro to compositional flow <a class="header-anchor" href="#Intro-to-compositional-flow" aria-label="Permalink to &quot;Intro to compositional flow {#Intro-to-compositional-flow}&quot;">​</a></h1><p>This is a simple conceptual example demonstrating how to solve compositional flow. This example uses a two-component water-CO2 system. Note that the default Peng-Robinson is not accurate for this system without adjustments to the parameters. However, the example demonstrates the conceptual workflow for getting started with compositional simulation.</p><h2 id="Set-up-mixture" tabindex="-1">Set up mixture <a class="header-anchor" href="#Set-up-mixture" aria-label="Permalink to &quot;Set up mixture {#Set-up-mixture}&quot;">​</a></h2><p>We load the external flash package and define a two-component H2O-CO2 system. The constructor for each species takes in molecular weight, critical pressure, critical temperature, critical volume, acentric factor given as strict SI. This means, for instance, that molar masses are given in kg/mole and not g/mole or kg/kmol.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> MultiComponentFlash</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">h2o </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> MolecularProperty</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.018015268</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">22.064e6</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">647.096</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5.595e-05</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.3442920843</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">co2 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> MolecularProperty</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0440098</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">7.3773e6</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">304.1282</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">9.412e-05</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.22394</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bic </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> zeros</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">mixture </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> MultiComponentMixture</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([h2o, co2], A_ij </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> bic, names </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;H2O&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;CO2&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">eos </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> GenericCubicEOS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(mixture, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">PengRobinson</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>MultiComponentFlash.GenericCubicEOS{MultiComponentFlash.PengRobinson, Float64, 2, Nothing}(MultiComponentFlash.PengRobinson(), MultiComponentFlash.MultiComponentMixture{Float64, 2}(&quot;UnnamedMixture&quot;, [&quot;H2O&quot;, &quot;CO2&quot;], (MultiComponentFlash.MolecularProperty{Float64}(0.018015268, 2.2064e7, 647.096, 5.595e-5, 0.3442920843), MultiComponentFlash.MolecularProperty{Float64}(0.0440098, 7.3773e6, 304.1282, 9.412e-5, 0.22394)), [0.0 0.0; 0.0 0.0]), 2.414213562373095, -0.41421356237309515, 0.457235529, 0.077796074, nothing)</span></span></code></pre></div><h2 id="Set-up-domain-and-wells" tabindex="-1">Set up domain and wells <a class="header-anchor" href="#Set-up-domain-and-wells" aria-label="Permalink to &quot;Set up domain and wells {#Set-up-domain-and-wells}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Jutul, JutulDarcy, GLMakie</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">nx </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 50</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ny </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">nz </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 20</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dims </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (nx, ny, nz)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">g </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> CartesianMesh</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(dims, (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">nc </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> number_of_cells</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(g)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Darcy, bar, kg, meter, Kelvin, day, sec </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> si_units</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:darcy</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:bar</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:kilogram</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:meter</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Kelvin</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:day</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:second</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">K </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> repeat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Darcy, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, nc)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">res </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reservoir_domain</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(g, porosity </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, permeability </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> K)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>DataDomain wrapping CartesianMesh (3D) with 50x1x20=1000 cells with 17 data fields added:</span></span>
<span class="line"><span>  1000 Cells</span></span>
<span class="line"><span>    :permeability =&gt; 3×1000 Matrix{Float64}</span></span>
<span class="line"><span>    :porosity =&gt; 1000 Vector{Float64}</span></span>
<span class="line"><span>    :rock_thermal_conductivity =&gt; 1000 Vector{Float64}</span></span>
<span class="line"><span>    :fluid_thermal_conductivity =&gt; 1000 Vector{Float64}</span></span>
<span class="line"><span>    :rock_density =&gt; 1000 Vector{Float64}</span></span>
<span class="line"><span>    :cell_centroids =&gt; 3×1000 Matrix{Float64}</span></span>
<span class="line"><span>    :volumes =&gt; 1000 Vector{Float64}</span></span>
<span class="line"><span>  1930 Faces</span></span>
<span class="line"><span>    :neighbors =&gt; 2×1930 Matrix{Int64}</span></span>
<span class="line"><span>    :areas =&gt; 1930 Vector{Float64}</span></span>
<span class="line"><span>    :normals =&gt; 3×1930 Matrix{Float64}</span></span>
<span class="line"><span>    :face_centroids =&gt; 3×1930 Matrix{Float64}</span></span>
<span class="line"><span>  3860 HalfFaces</span></span>
<span class="line"><span>    :half_face_cells =&gt; 3860 Vector{Int64}</span></span>
<span class="line"><span>    :half_face_faces =&gt; 3860 Vector{Int64}</span></span>
<span class="line"><span>  2140 BoundaryFaces</span></span>
<span class="line"><span>    :boundary_areas =&gt; 2140 Vector{Float64}</span></span>
<span class="line"><span>    :boundary_centroids =&gt; 3×2140 Matrix{Float64}</span></span>
<span class="line"><span>    :boundary_normals =&gt; 3×2140 Matrix{Float64}</span></span>
<span class="line"><span>    :boundary_neighbors =&gt; 2140 Vector{Int64}</span></span></code></pre></div><p>Set up a vertical well in the first corner, perforated in top layer</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">prod </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_well</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(g, K, [(nx, ny, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)], name </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Producer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>SimpleWell [Producer] (1 nodes, 0 segments, 1 perforations)</span></span></code></pre></div><p>Set up an injector in the opposite corner, perforated in bottom layer</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">inj </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_well</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(g, K, [(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, nz)], name </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Injector</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>SimpleWell [Injector] (1 nodes, 0 segments, 1 perforations)</span></span></code></pre></div><h2 id="Define-system-and-realize-on-grid" tabindex="-1">Define system and realize on grid <a class="header-anchor" href="#Define-system-and-realize-on-grid" aria-label="Permalink to &quot;Define system and realize on grid {#Define-system-and-realize-on-grid}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rhoLS </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 844.23</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">kg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">meter</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rhoVS </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 126.97</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">kg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">meter</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rhoS </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [rhoLS, rhoVS]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L, V </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> LiquidPhase</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">VaporPhase</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sys </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> MultiPhaseCompositionalSystemLV</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(eos, (L, V))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">model, parameters </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_reservoir_model</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(res, sys, wells </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [inj, prod]);</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">push!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">output_variables, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Saturations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">kr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> BrooksCoreyRelativePermeabilities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sys, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">model </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> replace_variables!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, RelativePermeabilities </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> kr)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">T0 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> fill</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">303.15</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Kelvin, nc)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">parameters[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Temperature</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> T0</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">state0 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_reservoir_state</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, Pressure </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 50</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar, OverallMoleFractions </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]);</span></span></code></pre></div><h2 id="Define-schedule" tabindex="-1">Define schedule <a class="header-anchor" href="#Define-schedule" aria-label="Permalink to &quot;Define schedule {#Define-schedule}&quot;">​</a></h2><p>5 year (5*365.24 days) simulation period</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dt0 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> fill</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">day, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">26</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dt1 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> fill</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10.0</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">day, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">180</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> cat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(dt0, dt1, dims </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rate_target </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> TotalRateTarget</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">9.5066e-06</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">meter</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sec)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">I_ctrl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> InjectorControl</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rate_target, [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], density </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> rhoVS)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bhp_target </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> BottomHolePressureTarget</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">50</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">P_ctrl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ProducerControl</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(bhp_target)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">controls </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Dict</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">controls[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Injector</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> I_ctrl</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">controls[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Producer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> P_ctrl</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">forces </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_reservoir_forces</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, control </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> controls)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ws, states </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(state0, model, dt, parameters </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> parameters, forces </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> forces);</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>\x1B[92;1mJutul:\x1B[0m Simulating 4 years, 52.15 weeks as 206 report steps</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Progress   1%|▍                                          |  ETA: 0:19:28\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 2/206 (0.11% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     6 iterations in 6.61 s (1.10 s each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress   5%|██▏                                        |  ETA: 0:03:47\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 10/206 (0.55% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     24 iterations in 6.71 s (279.77 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  11%|████▋                                      |  ETA: 0:01:38\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 22/206 (1.20% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     48 iterations in 6.82 s (142.08 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  15%|██████▋                                    |  ETA: 0:01:04\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 32/206 (4.71% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     79 iterations in 6.97 s (88.25 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  19%|████████▏                                  |  ETA: 0:00:51\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 39/206 (8.54% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     102 iterations in 7.10 s (69.63 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  22%|█████████▍                                 |  ETA: 0:00:43\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 45/206 (11.83% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     120 iterations in 7.22 s (60.13 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  24%|██████████▍                                |  ETA: 0:00:38\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 50/206 (14.57% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     135 iterations in 7.32 s (54.21 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  27%|███████████▍                               |  ETA: 0:00:34\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 55/206 (17.31% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     150 iterations in 7.43 s (49.53 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  29%|████████████▌                              |  ETA: 0:00:30\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 60/206 (20.04% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     167 iterations in 7.56 s (45.28 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  31%|█████████████▎                             |  ETA: 0:00:28\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 64/206 (22.23% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     180 iterations in 7.67 s (42.59 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  33%|██████████████▏                            |  ETA: 0:00:26\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 68/206 (24.42% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     194 iterations in 7.78 s (40.13 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  35%|███████████████                            |  ETA: 0:00:24\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 72/206 (26.62% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     207 iterations in 7.90 s (38.17 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  36%|███████████████▋                           |  ETA: 0:00:23\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 75/206 (28.26% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     219 iterations in 8.01 s (36.60 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  38%|████████████████▍                          |  ETA: 0:00:21\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 79/206 (30.45% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     232 iterations in 8.14 s (35.09 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  40%|█████████████████▎                         |  ETA: 0:00:20\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 83/206 (32.64% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     248 iterations in 8.30 s (33.46 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  42%|██████████████████▏                        |  ETA: 0:00:18\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 87/206 (34.83% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     259 iterations in 8.41 s (32.47 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  44%|██████████████████▉                        |  ETA: 0:00:17\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 91/206 (37.02% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     269 iterations in 8.52 s (31.67 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  46%|████████████████████                       |  ETA: 0:00:16\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 96/206 (39.76% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     279 iterations in 8.63 s (30.93 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  49%|█████████████████████                      |  ETA: 0:00:14\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 101/206 (42.50% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     290 iterations in 8.74 s (30.16 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  51%|██████████████████████                     |  ETA: 0:00:13\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 106/206 (45.24% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     302 iterations in 8.87 s (29.37 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  54%|███████████████████████                    |  ETA: 0:00:12\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 111/206 (47.97% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     312 iterations in 8.98 s (28.79 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  56%|███████████████████████▉                   |  ETA: 0:00:11\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 115/206 (50.16% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     320 iterations in 9.14 s (28.55 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  58%|████████████████████████▉                  |  ETA: 0:00:10\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 120/206 (52.90% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     331 iterations in 9.25 s (27.94 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  60%|██████████████████████████                 |  ETA: 0:00:09\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 125/206 (55.64% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     341 iterations in 9.35 s (27.43 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  63%|███████████████████████████                |  ETA: 0:00:08\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 130/206 (58.38% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     351 iterations in 9.46 s (26.96 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  65%|████████████████████████████               |  ETA: 0:00:08\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 135/206 (61.12% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     361 iterations in 9.57 s (26.51 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  68%|█████████████████████████████▏             |  ETA: 0:00:07\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 140/206 (63.86% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     371 iterations in 9.68 s (26.08 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  70%|██████████████████████████████▏            |  ETA: 0:00:06\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 145/206 (66.59% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     381 iterations in 9.78 s (25.66 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  72%|███████████████████████████████▏           |  ETA: 0:00:06\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 150/206 (69.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     391 iterations in 9.88 s (25.26 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  75%|████████████████████████████████▎          |  ETA: 0:00:05\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 155/206 (72.07% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     401 iterations in 9.98 s (24.89 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  77%|█████████████████████████████████▎         |  ETA: 0:00:04\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 160/206 (74.81% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     411 iterations in 10.09 s (24.55 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  80%|██████████████████████████████████▎        |  ETA: 0:00:04\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 165/206 (77.55% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     421 iterations in 10.19 s (24.21 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  82%|███████████████████████████████████▍       |  ETA: 0:00:03\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 170/206 (80.28% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     431 iterations in 10.29 s (23.89 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  85%|████████████████████████████████████▍      |  ETA: 0:00:03\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 175/206 (83.02% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     441 iterations in 10.40 s (23.58 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  87%|█████████████████████████████████████▍     |  ETA: 0:00:02\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 180/206 (85.76% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     451 iterations in 10.50 s (23.28 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  89%|██████████████████████████████████████▍    |  ETA: 0:00:02\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 185/206 (88.50% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     461 iterations in 10.60 s (23.00 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  92%|███████████████████████████████████████▌   |  ETA: 0:00:01\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 190/206 (91.24% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     471 iterations in 10.71 s (22.73 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  94%|████████████████████████████████████████▌  |  ETA: 0:00:01\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 195/206 (93.98% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     481 iterations in 10.81 s (22.47 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress  97%|█████████████████████████████████████████▊ |  ETA: 0:00:00\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 201/206 (97.26% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     493 iterations in 10.93 s (22.17 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Progress 100%|███████████████████████████████████████████| Time: 0:00:15\x1B[K</span></span>
<span class="line"><span>  Progress:  Solved step 206/206\x1B[K</span></span>
<span class="line"><span>  Stats:     505 iterations in 11.05 s (21.88 ms each)\x1B[K</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 206 steps │ 207 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   2.45146 │       2.43961 │  505 (0) │</span></span>
<span class="line"><span>│ Linearization  │   3.45631 │       3.43961 │  712 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   10.6019 │       10.5507 │ 2184 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   21.2039 │       21.1014 │ 4368 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬─────────┬────────────┬─────────╮</span></span>
<span class="line"><span>│ Timing type   │    Each │   Relative │   Total │</span></span>
<span class="line"><span>│               │      ms │ Percentage │       s │</span></span>
<span class="line"><span>├───────────────┼─────────┼────────────┼─────────┤</span></span>
<span class="line"><span>│ Properties    │  6.0533 │    27.67 % │  3.0569 │</span></span>
<span class="line"><span>│ Equations     │  3.0058 │    19.37 % │  2.1401 │</span></span>
<span class="line"><span>│ Assembly      │  1.4901 │     9.60 % │  1.0609 │</span></span>
<span class="line"><span>│ Linear solve  │  0.4534 │     2.07 % │  0.2290 │</span></span>
<span class="line"><span>│ Linear setup  │  1.1851 │     5.42 % │  0.5985 │</span></span>
<span class="line"><span>│ Precond apply │  0.0900 │     3.56 % │  0.3930 │</span></span>
<span class="line"><span>│ Update        │  0.8166 │     3.73 % │  0.4124 │</span></span>
<span class="line"><span>│ Convergence   │  1.9408 │    12.51 % │  1.3819 │</span></span>
<span class="line"><span>│ Input/Output  │  0.3542 │     0.66 % │  0.0733 │</span></span>
<span class="line"><span>│ Other         │  3.3714 │    15.41 % │  1.7025 │</span></span>
<span class="line"><span>├───────────────┼─────────┼────────────┼─────────┤</span></span>
<span class="line"><span>│ Total         │ 21.8784 │   100.00 % │ 11.0486 │</span></span>
<span class="line"><span>╰───────────────┴─────────┴────────────┴─────────╯</span></span></code></pre></div><h2 id="Once-the-simulation-is-done,-we-can-plot-the-states" tabindex="-1">Once the simulation is done, we can plot the states <a class="header-anchor" href="#Once-the-simulation-is-done,-we-can-plot-the-states" aria-label="Permalink to &quot;Once the simulation is done, we can plot the states {#Once-the-simulation-is-done,-we-can-plot-the-states}&quot;">​</a></h2><p>Note that this example is intended for static publication in the documentation. For interactive visualization you can use functions like <code>plot_interactive</code> to interactively visualize the states.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">z </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> states[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">end</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:OverallMoleFractions</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> plot_vertical</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x, t)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    data </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reshape</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x, (nx, nz))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    data </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data[:, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">end</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    fig, ax, plot </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> heatmap</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(data)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ax</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">title </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> t</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    Colorbar</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], plot)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    fig</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span></code></pre></div><h3 id="Plot-final-CO2-mole-fraction" tabindex="-1">Plot final CO2 mole fraction <a class="header-anchor" href="#Plot-final-CO2-mole-fraction" aria-label="Permalink to &quot;Plot final CO2 mole fraction {#Plot-final-CO2-mole-fraction}&quot;">​</a></h3><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_vertical</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(z, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;CO2&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="`+l+`" alt=""></p><h3 id="Plot-final-vapor-saturation" tabindex="-1">Plot final vapor saturation <a class="header-anchor" href="#Plot-final-vapor-saturation" aria-label="Permalink to &quot;Plot final vapor saturation {#Plot-final-vapor-saturation}&quot;">​</a></h3><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> states[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">end</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Saturations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_vertical</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sg, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Vapor saturation&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="`+e+`" alt=""></p><h3 id="Plot-final-pressure" tabindex="-1">Plot final pressure <a class="header-anchor" href="#Plot-final-pressure" aria-label="Permalink to &quot;Plot final pressure {#Plot-final-pressure}&quot;">​</a></h3><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> states[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">end</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Pressure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_vertical</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(p</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">./</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Pressure [bar]&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="`+t+'" alt=""></p><h3 id="Plot-in-interactive-viewer" tabindex="-1">Plot in interactive viewer <a class="header-anchor" href="#Plot-in-interactive-viewer" aria-label="Permalink to &quot;Plot in interactive viewer {#Plot-in-interactive-viewer}&quot;">​</a></h3><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, states, step </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> length</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(dt), key </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Saturations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="'+h+'" alt=""></p><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/co2_brine_2d_vertical.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/co2_brine_2d_vertical.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',40)]))}const y=a(k,[["render",r]]);export{b as __pageData,y as default};
