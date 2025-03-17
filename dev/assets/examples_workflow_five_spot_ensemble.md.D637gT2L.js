import{_ as e,c as i,o as p,az as n,j as a,a as l}from"./chunks/framework.yZYC2Rg6.js";const t="/JutulDarcy.jl/dev/assets/zjyzewu.DE8WqMTo.jpeg",h="/JutulDarcy.jl/dev/assets/oqqyilc.C1esbZIP.jpeg",k="/JutulDarcy.jl/dev/assets/tozihkp.B1cvO5k0.jpeg",r="/JutulDarcy.jl/dev/assets/twvysmy.Bb_G2Vi1.jpeg",v=JSON.parse('{"title":"Quarter-five-spot example","description":"","frontmatter":{},"headers":[],"relativePath":"examples/workflow/five_spot_ensemble.md","filePath":"examples/workflow/five_spot_ensemble.md","lastUpdated":null}'),E={name:"examples/workflow/five_spot_ensemble.md"},d={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},g={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.464ex"},xmlns:"http://www.w3.org/2000/svg",width:"5.42ex",height:"1.783ex",role:"img",focusable:"false",viewBox:"0 -583 2395.6 788","aria-hidden":"true"};function o(c,s,y,F,u,C){return p(),i("div",null,[s[4]||(s[4]=n(`<h1 id="Quarter-five-spot-example" tabindex="-1">Quarter-five-spot example <a class="header-anchor" href="#Quarter-five-spot-example" aria-label="Permalink to &quot;Quarter-five-spot example {#Quarter-five-spot-example}&quot;">​</a></h1><p>The quarter-five-spot is a standard test problem that simulates 1/4 of the five spot well pattern by assuming axial symmetry. The problem contains an injector in one corner and the producer in the opposing corner, with a significant volume of fluids injected into the domain.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy, Jutul</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">nx </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 50</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span></code></pre></div><h2 id="setup" tabindex="-1">Setup <a class="header-anchor" href="#setup" aria-label="Permalink to &quot;Setup&quot;">​</a></h2><p>We define a function that, for a given porosity field, computes a solution with an estimated permeability field. For assumptions and derivation of the specific form of the Kozeny-Carman relation used in this example, see <a href="https://doi.org/10.1017/9781108591416" target="_blank" rel="noreferrer">Lie, Knut-Andreas. An introduction to reservoir simulation using MATLAB/GNU Octave: User guide for the MATLAB Reservoir Simulation Toolbox (MRST). Cambridge University Press, 2019, Section 2.5.2</a></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> perm_kozeny_carman</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(Φ)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ((Φ</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.81</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">72</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Φ)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">);</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> simulate_qfs</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(porosity </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    Dx </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1000.0</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    Dz </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 10.0</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    Darcy </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 9.869232667160130e-13</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    Darcy, bar, kg, meter, Kelvin, day, sec </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> si_units</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:darcy</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:bar</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:kilogram</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:meter</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Kelvin</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:day</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:second</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    mesh </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> CartesianMesh</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">((nx, nx, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), (Dx, Dx, Dz))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    K </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> perm_kozeny_carman</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(porosity)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    domain </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reservoir_domain</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(mesh, permeability </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> K, porosity </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> porosity)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    Inj </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_vertical_well</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(domain, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, name </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Injector</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">);</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    Prod </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_vertical_well</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(domain, nx, nx, name </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Producer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">);</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    phases </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">LiquidPhase</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">VaporPhase</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    rhoLS </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1000.0</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">kg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">meter</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    rhoGS </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 700.0</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">kg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">meter</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    rhoS </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [rhoLS, rhoGS]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    sys </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ImmiscibleSystem</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(phases, reference_densities </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> rhoS)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    model, parameters </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_reservoir_model</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(domain, sys, wells </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [Inj, Prod])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    c </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-6</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-6</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ρ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ConstantCompressibilityDensities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(p_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 150</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar, density_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> rhoS, compressibility </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> c)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    kr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> BrooksCoreyRelativePermeabilities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sys, [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    replace_variables!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, PhaseMassDensities </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ρ, RelativePermeabilities </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> kr);</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    state0 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_reservoir_state</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, Pressure </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 150</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar, Saturations </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    dt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> repeat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">30.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">day, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">12</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    dt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> vcat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], dt)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    inj_rate </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Dx</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Dx</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Dz</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sum</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(dt) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># 1 PVI if average porosity is 0.2</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    rate_target </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> TotalRateTarget</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(inj_rate)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    I_ctrl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> InjectorControl</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rate_target, [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], density </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> rhoGS)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    bhp_target </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> BottomHolePressureTarget</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">50</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    P_ctrl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ProducerControl</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(bhp_target)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    controls </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Dict</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    controls[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Injector</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> I_ctrl</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    controls[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Producer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> P_ctrl</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    forces </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_reservoir_forces</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, control </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> controls)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(state0, model, dt, parameters </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> parameters, forces </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> forces)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>simulate_qfs (generic function with 2 methods)</span></span></code></pre></div><h2 id="Simulate-base-case" tabindex="-1">Simulate base case <a class="header-anchor" href="#Simulate-base-case" aria-label="Permalink to &quot;Simulate base case {#Simulate-base-case}&quot;">​</a></h2><p>This will give the solution with uniform porosity of 0.2.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ws, states, report_time </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_qfs</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">();</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 141 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │    3.3252 │       2.90071 │  409 (0) │</span></span>
<span class="line"><span>│ Linearization  │   4.47154 │       3.90071 │  550 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   10.2358 │       8.92908 │ 1259 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   20.4715 │       17.8582 │ 2518 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2447 │     5.13 % │ 0.1001 │</span></span>
<span class="line"><span>│ Equations     │ 0.2381 │     6.71 % │ 0.1310 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2453 │     6.91 % │ 0.1349 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3306 │     6.92 % │ 0.1352 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.0508 │    42.95 % │ 0.8388 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1908 │    24.60 % │ 0.4804 │</span></span>
<span class="line"><span>│ Update        │ 0.0911 │     1.91 % │ 0.0373 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0881 │     2.48 % │ 0.0485 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0522 │     0.38 % │ 0.0074 │</span></span>
<span class="line"><span>│ Other         │ 0.0968 │     2.03 % │ 0.0396 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.7752 │   100.00 % │ 1.9530 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span></code></pre></div><h3 id="Plot-the-solution-of-the-base-case" tabindex="-1">Plot the solution of the base case <a class="header-anchor" href="#Plot-the-solution-of-the-base-case" aria-label="Permalink to &quot;Plot the solution of the base case {#Plot-the-solution-of-the-base-case}&quot;">​</a></h3>`,12)),a("p",null,[s[2]||(s[2]=l("We observe a radial flow pattern initially, before coning occurs near the producer well once the fluid has reached the opposite corner. The uniform permeability and porosity gives axial symmetry at ")),a("mjx-container",d,[(p(),i("svg",g,s[0]||(s[0]=[n('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="mi"><path data-c="1D465" d="M52 289Q59 331 106 386T222 442Q257 442 286 424T329 379Q371 442 430 442Q467 442 494 420T522 361Q522 332 508 314T481 292T458 288Q439 288 427 299T415 328Q415 374 465 391Q454 404 425 404Q412 404 406 402Q368 386 350 336Q290 115 290 78Q290 50 306 38T341 26Q378 26 414 59T463 140Q466 150 469 151T485 153H489Q504 153 504 145Q504 144 502 134Q486 77 440 33T333 -11Q263 -11 227 52Q186 -10 133 -10H127Q78 -10 57 16T35 71Q35 103 54 123T99 143Q142 143 142 101Q142 81 130 66T107 46T94 41L91 40Q91 39 97 36T113 29T132 26Q168 26 194 71Q203 87 217 139T245 247T261 313Q266 340 266 352Q266 380 251 392T217 404Q177 404 142 372T93 290Q91 281 88 280T72 278H58Q52 284 52 289Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(849.8,0)"><path data-c="3D" d="M56 347Q56 360 70 367H707Q722 359 722 347Q722 336 708 328L390 327H72Q56 332 56 347ZM56 153Q56 168 72 173H708Q722 163 722 153Q722 140 707 133H70Q56 140 56 153Z" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(1905.6,0)"><path data-c="1D466" d="M21 287Q21 301 36 335T84 406T158 442Q199 442 224 419T250 355Q248 336 247 334Q247 331 231 288T198 191T182 105Q182 62 196 45T238 27Q261 27 281 38T312 61T339 94Q339 95 344 114T358 173T377 247Q415 397 419 404Q432 431 462 431Q475 431 483 424T494 412T496 403Q496 390 447 193T391 -23Q363 -106 294 -155T156 -205Q111 -205 77 -183T43 -117Q43 -95 50 -80T69 -58T89 -48T106 -45Q150 -45 150 -87Q150 -107 138 -122T115 -142T102 -147L99 -148Q101 -153 118 -160T152 -167H160Q177 -167 186 -165Q219 -156 247 -127T290 -65T313 -9T321 21L315 17Q309 13 296 6T270 -6Q250 -11 231 -11Q185 -11 150 11T104 82Q103 89 103 113Q103 170 138 262T173 379Q173 380 173 381Q173 390 173 393T169 400T158 404H154Q131 404 112 385T82 344T65 302T57 280Q55 278 41 278H27Q21 284 21 287Z" style="stroke-width:3;"></path></g></g></g>',1)]))),s[1]||(s[1]=a("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[a("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[a("mi",null,"x"),a("mo",null,"="),a("mi",null,"y")])],-1))]),s[3]||(s[3]=l("."))]),s[5]||(s[5]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GLMakie</span></span>
<span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">to_2d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reshape</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">vec</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x), nx, nx)</span></span>
<span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">get_sat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(state) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> to_2d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(state[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Saturations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">nt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> length</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(report_time)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">h </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> nothing</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">h </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> contourf!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_sat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(states[nt</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">÷</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">h </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> contourf!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_sat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(states[nt]))</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Colorbar</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">end</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], h)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+t+`" alt=""></p><h2 id="Create-10-realizations" tabindex="-1">Create 10 realizations <a class="header-anchor" href="#Create-10-realizations" aria-label="Permalink to &quot;Create 10 realizations {#Create-10-realizations}&quot;">​</a></h2><p>We create a small set of realizations of the same model, with porosity that is uniformly varying between 0.05 and 0.3. This is not especially sophisticated geostatistics - for a more realistic approach, take a look at <a href="https://juliaearth.github.io/GeoStats.jl" target="_blank" rel="noreferrer">GeoStats.jl</a>. The main idea is to get significantly different flow patterns as the porosity and permeability changes.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">N </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 10</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">saturations </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> []</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">wells </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> []</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">report_step </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> nt</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">N</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    poro </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.05</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> .+</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.25</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">rand</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(Float64, (nx</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">nx))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ws_i, states_i, rt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_qfs</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(poro)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    push!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(wells, ws_i)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    push!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(saturations, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_sat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(states_i[report_step]))</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬───────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │     Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 158 ministeps │  (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼───────────┤</span></span>
<span class="line"><span>│ Newton         │   4.45528 │       3.46835 │  548 (15) │</span></span>
<span class="line"><span>│ Linearization  │   5.73984 │       4.46835 │  706 (16) │</span></span>
<span class="line"><span>│ Linear solver  │   11.2927 │       8.79114 │ 1389 (15) │</span></span>
<span class="line"><span>│ Precond apply  │   22.5854 │       17.5823 │ 2778 (30) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴───────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2461 │     4.48 % │ 0.1349 │</span></span>
<span class="line"><span>│ Equations     │ 0.2580 │     6.05 % │ 0.1821 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2505 │     5.87 % │ 0.1769 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2910 │     5.30 % │ 0.1595 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2002 │    40.04 % │ 1.2057 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1869 │    17.25 % │ 0.5193 │</span></span>
<span class="line"><span>│ Update        │ 0.0941 │     1.71 % │ 0.0516 │</span></span>
<span class="line"><span>│ Convergence   │ 0.6661 │    15.62 % │ 0.4702 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0508 │     0.27 % │ 0.0080 │</span></span>
<span class="line"><span>│ Other         │ 0.1874 │     3.41 % │ 0.1027 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.4943 │   100.00 % │ 3.0109 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 152 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   4.02439 │       3.25658 │  495 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.26016 │       4.25658 │  647 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   18.9593 │       15.3421 │ 2332 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   37.9187 │       30.6842 │ 4664 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2516 │     4.44 % │ 0.1246 │</span></span>
<span class="line"><span>│ Equations     │ 0.2713 │     6.26 % │ 0.1755 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2593 │     5.99 % │ 0.1678 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.4374 │     7.72 % │ 0.2165 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2427 │    39.60 % │ 1.1101 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1800 │    29.95 % │ 0.8397 │</span></span>
<span class="line"><span>│ Update        │ 0.1037 │     1.83 % │ 0.0513 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0963 │     2.22 % │ 0.0623 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0576 │     0.31 % │ 0.0087 │</span></span>
<span class="line"><span>│ Other         │ 0.0944 │     1.67 % │ 0.0467 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.6634 │   100.00 % │ 2.8034 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 154 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   4.06504 │       3.24675 │  500 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.31707 │       4.24675 │  654 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   18.5285 │       14.7987 │ 2279 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   37.0569 │       29.5974 │ 4558 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2486 │     4.42 % │ 0.1243 │</span></span>
<span class="line"><span>│ Equations     │ 0.2643 │     6.15 % │ 0.1728 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2522 │     5.87 % │ 0.1649 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.4334 │     7.71 % │ 0.2167 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2911 │    40.76 % │ 1.1455 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1801 │    29.21 % │ 0.8207 │</span></span>
<span class="line"><span>│ Update        │ 0.0969 │     1.72 % │ 0.0485 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0948 │     2.21 % │ 0.0620 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0546 │     0.30 % │ 0.0084 │</span></span>
<span class="line"><span>│ Other         │ 0.0925 │     1.65 % │ 0.0462 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.6204 │   100.00 % │ 2.8102 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 153 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   4.21138 │       3.38562 │  518 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.45528 │       4.38562 │  671 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   19.4472 │        15.634 │ 2392 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   38.8943 │        31.268 │ 4784 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2542 │     4.40 % │ 0.1317 │</span></span>
<span class="line"><span>│ Equations     │ 0.2883 │     6.47 % │ 0.1934 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2642 │     5.93 % │ 0.1773 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.4582 │     7.93 % │ 0.2374 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.3186 │    40.15 % │ 1.2010 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1807 │    28.89 % │ 0.8643 │</span></span>
<span class="line"><span>│ Update        │ 0.1069 │     1.85 % │ 0.0554 │</span></span>
<span class="line"><span>│ Convergence   │ 0.1063 │     2.38 % │ 0.0713 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0565 │     0.29 % │ 0.0086 │</span></span>
<span class="line"><span>│ Other         │ 0.0986 │     1.71 % │ 0.0511 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.7750 │   100.00 % │ 2.9914 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬────────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │      Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 213 ministeps │   (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼────────────┤</span></span>
<span class="line"><span>│ Newton         │    8.6748 │       5.00939 │ 1067 (330) │</span></span>
<span class="line"><span>│ Linearization  │   10.4065 │       6.00939 │ 1280 (352) │</span></span>
<span class="line"><span>│ Linear solver  │   16.3577 │       9.44601 │ 2012 (356) │</span></span>
<span class="line"><span>│ Precond apply  │   32.7154 │        18.892 │ 4024 (712) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴────────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2418 │     6.07 % │ 0.2580 │</span></span>
<span class="line"><span>│ Equations     │ 0.2437 │     7.34 % │ 0.3120 │</span></span>
<span class="line"><span>│ Assembly      │ 0.3096 │     9.32 % │ 0.3963 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2556 │     6.42 % │ 0.2727 │</span></span>
<span class="line"><span>│ Linear setup  │ 1.8676 │    46.88 % │ 1.9928 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1844 │    17.46 % │ 0.7422 │</span></span>
<span class="line"><span>│ Update        │ 0.0866 │     2.18 % │ 0.0924 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0843 │     2.54 % │ 0.1079 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0417 │     0.21 % │ 0.0089 │</span></span>
<span class="line"><span>│ Other         │ 0.0630 │     1.58 % │ 0.0672 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 3.9834 │   100.00 % │ 4.2503 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 148 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │    3.6748 │       3.05405 │  452 (0) │</span></span>
<span class="line"><span>│ Linearization  │   4.87805 │       4.05405 │  600 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   15.3902 │       12.7905 │ 1893 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   30.7805 │       25.5811 │ 3786 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2508 │     4.53 % │ 0.1133 │</span></span>
<span class="line"><span>│ Equations     │ 0.2694 │     6.46 % │ 0.1616 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2540 │     6.09 % │ 0.1524 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.4184 │     7.56 % │ 0.1891 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.3067 │    41.68 % │ 1.0426 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1813 │    27.43 % │ 0.6863 │</span></span>
<span class="line"><span>│ Update        │ 0.1014 │     1.83 % │ 0.0458 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0962 │     2.31 % │ 0.0577 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0551 │     0.33 % │ 0.0082 │</span></span>
<span class="line"><span>│ Other         │ 0.0988 │     1.79 % │ 0.0447 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.5350 │   100.00 % │ 2.5018 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬─────────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │       Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 322 ministeps │    (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼─────────────┤</span></span>
<span class="line"><span>│ Newton         │   16.6829 │       6.37267 │  2052 (930) │</span></span>
<span class="line"><span>│ Linearization  │   19.3008 │       7.37267 │  2374 (992) │</span></span>
<span class="line"><span>│ Linear solver  │   24.9756 │       9.54037 │  3072 (963) │</span></span>
<span class="line"><span>│ Precond apply  │   49.9512 │       19.0807 │ 6144 (1926) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴─────────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2391 │     6.43 % │ 0.4906 │</span></span>
<span class="line"><span>│ Equations     │ 0.2671 │     8.31 % │ 0.6342 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2407 │     7.49 % │ 0.5715 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2364 │     6.36 % │ 0.4852 │</span></span>
<span class="line"><span>│ Linear setup  │ 1.8579 │    49.99 % │ 3.8125 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1861 │    14.99 % │ 1.1437 │</span></span>
<span class="line"><span>│ Update        │ 0.0829 │     2.23 % │ 0.1701 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0829 │     2.58 % │ 0.1967 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0352 │     0.15 % │ 0.0113 │</span></span>
<span class="line"><span>│ Other         │ 0.0543 │     1.46 % │ 0.1114 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 3.7169 │   100.00 % │ 7.6271 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 158 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   4.21138 │       3.27848 │  518 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.49593 │       4.27848 │  676 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   21.5203 │       16.7532 │ 2647 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   43.0407 │       33.5063 │ 5294 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2483 │     4.25 % │ 0.1286 │</span></span>
<span class="line"><span>│ Equations     │ 0.2630 │     5.87 % │ 0.1778 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2520 │     5.63 % │ 0.1704 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.4636 │     7.93 % │ 0.2402 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.3019 │    39.38 % │ 1.1924 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1794 │    31.36 % │ 0.9496 │</span></span>
<span class="line"><span>│ Update        │ 0.0959 │     1.64 % │ 0.0497 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0949 │     2.12 % │ 0.0641 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0502 │     0.26 % │ 0.0079 │</span></span>
<span class="line"><span>│ Other         │ 0.0907 │     1.55 % │ 0.0470 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.8450 │   100.00 % │ 3.0277 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬────────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │      Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 173 ministeps │   (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼────────────┤</span></span>
<span class="line"><span>│ Newton         │   5.34959 │       3.80347 │   658 (15) │</span></span>
<span class="line"><span>│ Linearization  │    6.7561 │       4.80347 │   831 (16) │</span></span>
<span class="line"><span>│ Linear solver  │   26.0894 │       18.5491 │  3209 (91) │</span></span>
<span class="line"><span>│ Precond apply  │   52.1789 │       37.0983 │ 6418 (182) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴────────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2502 │     4.27 % │ 0.1647 │</span></span>
<span class="line"><span>│ Equations     │ 0.2723 │     5.87 % │ 0.2263 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2567 │     5.53 % │ 0.2133 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.4628 │     7.89 % │ 0.3045 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.3289 │    39.72 % │ 1.5324 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1804 │    30.01 % │ 1.1577 │</span></span>
<span class="line"><span>│ Update        │ 0.1019 │     1.74 % │ 0.0670 │</span></span>
<span class="line"><span>│ Convergence   │ 0.1014 │     2.19 % │ 0.0843 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.2968 │     1.33 % │ 0.0514 │</span></span>
<span class="line"><span>│ Other         │ 0.0861 │     1.47 % │ 0.0567 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.8636 │   100.00 % │ 3.8583 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬─────────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │       Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 362 ministeps │    (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼─────────────┤</span></span>
<span class="line"><span>│ Newton         │    20.065 │       6.81768 │ 2468 (1185) │</span></span>
<span class="line"><span>│ Linearization  │   23.0081 │       7.81768 │ 2830 (1264) │</span></span>
<span class="line"><span>│ Linear solver  │    30.187 │       10.2569 │ 3713 (1218) │</span></span>
<span class="line"><span>│ Precond apply  │    60.374 │       20.5138 │ 7426 (2436) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴─────────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2377 │     6.52 % │ 0.5867 │</span></span>
<span class="line"><span>│ Equations     │ 0.2452 │     7.72 % │ 0.6940 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2674 │     8.41 % │ 0.7566 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2338 │     6.42 % │ 0.5770 │</span></span>
<span class="line"><span>│ Linear setup  │ 1.7967 │    49.31 % │ 4.4342 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1840 │    15.20 % │ 1.3666 │</span></span>
<span class="line"><span>│ Update        │ 0.0819 │     2.25 % │ 0.2022 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0831 │     2.61 % │ 0.2351 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0347 │     0.14 % │ 0.0126 │</span></span>
<span class="line"><span>│ Other         │ 0.0517 │     1.42 % │ 0.1275 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 3.6436 │   100.00 % │ 8.9924 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span></code></pre></div><h3 id="Plot-the-oil-rate-at-the-producer-over-the-ensemble" tabindex="-1">Plot the oil rate at the producer over the ensemble <a class="header-anchor" href="#Plot-the-oil-rate-at-the-producer-over-the-ensemble" aria-label="Permalink to &quot;Plot the oil rate at the producer over the ensemble {#Plot-the-oil-rate-at-the-producer-over-the-ensemble}&quot;">​</a></h3><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Statistics</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">N</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ws </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> wells[i]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    q </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ws[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Producer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:orat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, report_time, q)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">xlims!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">mean</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(report_time), report_time[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">end</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]])</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ylims!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0075</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+h+`" alt=""></p><h3 id="Plot-the-average-saturation-over-the-ensemble" tabindex="-1">Plot the average saturation over the ensemble <a class="header-anchor" href="#Plot-the-average-saturation-over-the-ensemble" aria-label="Permalink to &quot;Plot the average saturation over the ensemble {#Plot-the-average-saturation-over-the-ensemble}&quot;">​</a></h3><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">avg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> mean</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(saturations)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">h </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> nothing</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">h </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> contourf!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, avg)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+k+`" alt=""></p><h3 id="Plot-the-isocontour-lines-over-the-ensemble" tabindex="-1">Plot the isocontour lines over the ensemble <a class="header-anchor" href="#Plot-the-isocontour-lines-over-the-ensemble" aria-label="Permalink to &quot;Plot the isocontour lines over the ensemble {#Plot-the-isocontour-lines-over-the-ensemble}&quot;">​</a></h3><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">h </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> nothing</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> s </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> saturations</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    contour!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, s, levels </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+r+'" alt=""></p><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/workflow/five_spot_ensemble.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/workflow/five_spot_ensemble.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>This example took 52.176843942 seconds to complete.</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',20))])}const A=e(E,[["render",o]]);export{v as __pageData,A as default};
