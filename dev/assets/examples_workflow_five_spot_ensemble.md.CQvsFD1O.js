import{_ as e,c as i,o as p,az as n,j as a,a as l}from"./chunks/framework.yZYC2Rg6.js";const t="/JutulDarcy.jl/dev/assets/zjyzewu.DE8WqMTo.jpeg",h="/JutulDarcy.jl/dev/assets/oqqyilc.Cr3NvpY3.jpeg",k="/JutulDarcy.jl/dev/assets/tozihkp.Bu-Otzy-.jpeg",r="/JutulDarcy.jl/dev/assets/twvysmy.Cw5a7PnF.jpeg",v=JSON.parse('{"title":"Quarter-five-spot example","description":"","frontmatter":{},"headers":[],"relativePath":"examples/workflow/five_spot_ensemble.md","filePath":"examples/workflow/five_spot_ensemble.md","lastUpdated":null}'),E={name:"examples/workflow/five_spot_ensemble.md"},d={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},g={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.464ex"},xmlns:"http://www.w3.org/2000/svg",width:"5.42ex",height:"1.783ex",role:"img",focusable:"false",viewBox:"0 -583 2395.6 788","aria-hidden":"true"};function o(c,s,y,F,u,C){return p(),i("div",null,[s[4]||(s[4]=n(`<h1 id="Quarter-five-spot-example" tabindex="-1">Quarter-five-spot example <a class="header-anchor" href="#Quarter-five-spot-example" aria-label="Permalink to &quot;Quarter-five-spot example {#Quarter-five-spot-example}&quot;">​</a></h1><p>The quarter-five-spot is a standard test problem that simulates 1/4 of the five spot well pattern by assuming axial symmetry. The problem contains an injector in one corner and the producer in the opposing corner, with a significant volume of fluids injected into the domain.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy, Jutul</span></span>
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
<span class="line"><span>│ Properties    │ 0.2383 │     5.47 % │ 0.0974 │</span></span>
<span class="line"><span>│ Equations     │ 0.2076 │     6.41 % │ 0.1142 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2548 │     7.87 % │ 0.1401 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2880 │     6.62 % │ 0.1178 │</span></span>
<span class="line"><span>│ Linear setup  │ 1.7969 │    41.29 % │ 0.7349 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1889 │    26.73 % │ 0.4757 │</span></span>
<span class="line"><span>│ Update        │ 0.0741 │     1.70 % │ 0.0303 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0677 │     2.09 % │ 0.0372 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0340 │     0.27 % │ 0.0048 │</span></span>
<span class="line"><span>│ Other         │ 0.0670 │     1.54 % │ 0.0274 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.3519 │   100.00 % │ 1.7799 │</span></span>
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
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 153 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │       4.0 │       3.21569 │  492 (0) │</span></span>
<span class="line"><span>│ Linearization  │    5.2439 │       4.21569 │  645 (0) │</span></span>
<span class="line"><span>│ Linear solver  │      12.0 │       9.64706 │ 1476 (0) │</span></span>
<span class="line"><span>│ Precond apply  │      24.0 │       19.2941 │ 2952 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2399 │     5.07 % │ 0.1181 │</span></span>
<span class="line"><span>│ Equations     │ 0.2365 │     6.55 % │ 0.1526 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2634 │     7.30 % │ 0.1699 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2834 │     5.99 % │ 0.1394 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2122 │    46.75 % │ 1.0884 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1821 │    23.09 % │ 0.5377 │</span></span>
<span class="line"><span>│ Update        │ 0.0800 │     1.69 % │ 0.0394 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0715 │     1.98 % │ 0.0461 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0330 │     0.22 % │ 0.0050 │</span></span>
<span class="line"><span>│ Other         │ 0.0644 │     1.36 % │ 0.0317 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.7321 │   100.00 % │ 2.3282 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 146 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   3.66667 │       3.08904 │  451 (0) │</span></span>
<span class="line"><span>│ Linearization  │   4.85366 │       4.08904 │  597 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   15.8862 │       13.3836 │ 1954 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   31.7724 │       26.7671 │ 3908 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2417 │     4.58 % │ 0.1090 │</span></span>
<span class="line"><span>│ Equations     │ 0.2373 │     5.95 % │ 0.1416 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2653 │     6.66 % │ 0.1584 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3597 │     6.82 % │ 0.1622 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2051 │    41.80 % │ 0.9945 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1784 │    29.30 % │ 0.6973 │</span></span>
<span class="line"><span>│ Update        │ 0.0833 │     1.58 % │ 0.0376 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0727 │     1.82 % │ 0.0434 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0354 │     0.22 % │ 0.0052 │</span></span>
<span class="line"><span>│ Other         │ 0.0670 │     1.27 % │ 0.0302 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.2758 │   100.00 % │ 2.3794 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 149 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   3.95935 │       3.26846 │  487 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.17073 │       4.26846 │  636 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   18.6341 │       15.3826 │ 2292 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   37.2683 │       30.7651 │ 4584 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2397 │     4.24 % │ 0.1167 │</span></span>
<span class="line"><span>│ Equations     │ 0.4081 │     9.42 % │ 0.2595 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2571 │     5.93 % │ 0.1635 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3829 │     6.77 % │ 0.1865 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2431 │    39.64 % │ 1.0924 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1790 │    29.78 % │ 0.8207 │</span></span>
<span class="line"><span>│ Update        │ 0.0758 │     1.34 % │ 0.0369 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0693 │     1.60 % │ 0.0440 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0348 │     0.19 % │ 0.0052 │</span></span>
<span class="line"><span>│ Other         │ 0.0621 │     1.10 % │ 0.0302 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.6585 │   100.00 % │ 2.7557 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬────────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │      Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 174 ministeps │   (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼────────────┤</span></span>
<span class="line"><span>│ Newton         │   5.81301 │        4.1092 │   715 (60) │</span></span>
<span class="line"><span>│ Linearization  │   7.22764 │        5.1092 │   889 (64) │</span></span>
<span class="line"><span>│ Linear solver  │   13.1382 │       9.28736 │  1616 (60) │</span></span>
<span class="line"><span>│ Precond apply  │   26.2764 │       18.5747 │ 3232 (120) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴────────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2375 │     4.90 % │ 0.1698 │</span></span>
<span class="line"><span>│ Equations     │ 0.2249 │     5.77 % │ 0.1999 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2565 │     6.58 % │ 0.2280 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2413 │     4.98 % │ 0.1725 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.0890 │    43.10 % │ 1.4936 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1842 │    17.17 % │ 0.5952 │</span></span>
<span class="line"><span>│ Update        │ 0.0743 │     1.53 % │ 0.0531 │</span></span>
<span class="line"><span>│ Convergence   │ 0.5162 │    13.24 % │ 0.4589 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0309 │     0.15 % │ 0.0054 │</span></span>
<span class="line"><span>│ Other         │ 0.1250 │     2.58 % │ 0.0894 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.8473 │   100.00 % │ 3.4658 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 156 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   4.06504 │       3.20513 │  500 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.33333 │       4.20513 │  656 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   17.5935 │       13.8718 │ 2164 (0) │</span></span>
<span class="line"><span>│ Precond apply  │    35.187 │       27.7436 │ 4328 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2433 │     4.45 % │ 0.1217 │</span></span>
<span class="line"><span>│ Equations     │ 0.2371 │     5.69 % │ 0.1556 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2630 │     6.31 % │ 0.1726 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3604 │     6.59 % │ 0.1802 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2293 │    40.75 % │ 1.1146 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1792 │    28.35 % │ 0.7755 │</span></span>
<span class="line"><span>│ Update        │ 0.2606 │     4.76 % │ 0.1303 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0717 │     1.72 % │ 0.0470 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0343 │     0.20 % │ 0.0053 │</span></span>
<span class="line"><span>│ Other         │ 0.0648 │     1.18 % │ 0.0324 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.4703 │   100.00 % │ 2.7352 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬────────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │      Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 188 ministeps │   (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼────────────┤</span></span>
<span class="line"><span>│ Newton         │   6.74797 │       4.41489 │  830 (105) │</span></span>
<span class="line"><span>│ Linearization  │   8.27642 │       5.41489 │ 1018 (112) │</span></span>
<span class="line"><span>│ Linear solver  │   16.0813 │       10.5213 │ 1978 (149) │</span></span>
<span class="line"><span>│ Precond apply  │   32.1626 │       21.0426 │ 3956 (298) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴────────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2359 │     5.58 % │ 0.1958 │</span></span>
<span class="line"><span>│ Equations     │ 0.2210 │     6.41 % │ 0.2250 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2546 │     7.38 % │ 0.2592 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2494 │     5.90 % │ 0.2070 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.0818 │    49.21 % │ 1.7279 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1818 │    20.49 % │ 0.7194 │</span></span>
<span class="line"><span>│ Update        │ 0.0731 │     1.73 % │ 0.0606 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0694 │     2.01 % │ 0.0706 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0294 │     0.16 % │ 0.0055 │</span></span>
<span class="line"><span>│ Other         │ 0.0489 │     1.16 % │ 0.0406 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.2308 │   100.00 % │ 3.5116 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 148 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   3.88618 │       3.22973 │  478 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.08943 │       4.22973 │  626 (0) │</span></span>
<span class="line"><span>│ Linear solver  │    18.439 │       15.3243 │ 2268 (0) │</span></span>
<span class="line"><span>│ Precond apply  │    36.878 │       30.6486 │ 4536 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2402 │     4.31 % │ 0.1148 │</span></span>
<span class="line"><span>│ Equations     │ 0.2302 │     5.41 % │ 0.1441 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2601 │     6.11 % │ 0.1628 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3843 │     6.89 % │ 0.1837 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.3700 │    42.51 % │ 1.1329 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1785 │    30.38 % │ 0.8096 │</span></span>
<span class="line"><span>│ Update        │ 0.0788 │     1.41 % │ 0.0376 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0701 │     1.65 % │ 0.0439 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0345 │     0.19 % │ 0.0051 │</span></span>
<span class="line"><span>│ Other         │ 0.0636 │     1.14 % │ 0.0304 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.5753 │   100.00 % │ 2.6650 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬────────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │      Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 184 ministeps │   (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼────────────┤</span></span>
<span class="line"><span>│ Newton         │   6.30894 │       4.21739 │  776 (105) │</span></span>
<span class="line"><span>│ Linearization  │   7.80488 │       5.21739 │  960 (112) │</span></span>
<span class="line"><span>│ Linear solver  │   13.3659 │       8.93478 │ 1644 (105) │</span></span>
<span class="line"><span>│ Precond apply  │   26.7317 │       17.8696 │ 3288 (210) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴────────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2365 │     5.79 % │ 0.1835 │</span></span>
<span class="line"><span>│ Equations     │ 0.2232 │     6.76 % │ 0.2142 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2553 │     7.73 % │ 0.2450 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2364 │     5.79 % │ 0.1835 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.0271 │    49.63 % │ 1.5730 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1834 │    19.03 % │ 0.6031 │</span></span>
<span class="line"><span>│ Update        │ 0.0736 │     1.80 % │ 0.0571 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0687 │     2.08 % │ 0.0660 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0292 │     0.17 % │ 0.0054 │</span></span>
<span class="line"><span>│ Other         │ 0.0501 │     1.23 % │ 0.0388 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.0847 │   100.00 % │ 3.1697 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 155 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   4.12195 │       3.27097 │  507 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.38211 │       4.27097 │  662 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   15.4309 │       12.2452 │ 1898 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   30.8618 │       24.4903 │ 3796 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2421 │     4.70 % │ 0.1227 │</span></span>
<span class="line"><span>│ Equations     │ 0.2259 │     5.72 % │ 0.1496 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2544 │     6.44 % │ 0.1684 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.4009 │     7.78 % │ 0.2033 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2706 │    44.06 % │ 1.1512 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1835 │    26.67 % │ 0.6967 │</span></span>
<span class="line"><span>│ Update        │ 0.0752 │     1.46 % │ 0.0381 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0693 │     1.76 % │ 0.0459 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0340 │     0.20 % │ 0.0053 │</span></span>
<span class="line"><span>│ Other         │ 0.0625 │     1.21 % │ 0.0317 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.1535 │   100.00 % │ 2.6128 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 152 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   4.11382 │       3.32895 │  506 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.34959 │       4.32895 │  658 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   12.2602 │       9.92105 │ 1508 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   24.5203 │       19.8421 │ 3016 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2369 │     5.00 % │ 0.1199 │</span></span>
<span class="line"><span>│ Equations     │ 0.2225 │     6.10 % │ 0.1464 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2558 │     7.02 % │ 0.1683 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2805 │     5.92 % │ 0.1419 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2129 │    46.69 % │ 1.1197 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1823 │    22.92 % │ 0.5497 │</span></span>
<span class="line"><span>│ Update        │ 0.0747 │     1.58 % │ 0.0378 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0685 │     1.88 % │ 0.0451 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.2577 │     1.63 % │ 0.0392 │</span></span>
<span class="line"><span>│ Other         │ 0.0599 │     1.26 % │ 0.0303 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.7397 │   100.00 % │ 2.3983 │</span></span>
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
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+r+'" alt=""></p><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/workflow/five_spot_ensemble.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/workflow/five_spot_ensemble.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>This example took 37.524909499 seconds to complete.</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',20))])}const A=e(E,[["render",o]]);export{v as __pageData,A as default};
