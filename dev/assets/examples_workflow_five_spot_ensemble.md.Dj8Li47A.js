import{_ as e,c as i,o as p,az as n,j as a,a as l}from"./chunks/framework.yZYC2Rg6.js";const t="/JutulDarcy.jl/dev/assets/zjyzewu.0Tw0a6BN.jpeg",h="/JutulDarcy.jl/dev/assets/oqqyilc.DSWwPn7K.jpeg",k="/JutulDarcy.jl/dev/assets/tozihkp.DBZyq5Ps.jpeg",r="/JutulDarcy.jl/dev/assets/twvysmy.MM6aFfgc.jpeg",v=JSON.parse('{"title":"Quarter-five-spot example","description":"","frontmatter":{},"headers":[],"relativePath":"examples/workflow/five_spot_ensemble.md","filePath":"examples/workflow/five_spot_ensemble.md","lastUpdated":null}'),E={name:"examples/workflow/five_spot_ensemble.md"},d={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},g={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.464ex"},xmlns:"http://www.w3.org/2000/svg",width:"5.42ex",height:"1.783ex",role:"img",focusable:"false",viewBox:"0 -583 2395.6 788","aria-hidden":"true"};function o(c,s,y,F,u,C){return p(),i("div",null,[s[4]||(s[4]=n(`<h1 id="Quarter-five-spot-example" tabindex="-1">Quarter-five-spot example <a class="header-anchor" href="#Quarter-five-spot-example" aria-label="Permalink to &quot;Quarter-five-spot example {#Quarter-five-spot-example}&quot;">​</a></h1><p>The quarter-five-spot is a standard test problem that simulates 1/4 of the five spot well pattern by assuming axial symmetry. The problem contains an injector in one corner and the producer in the opposing corner, with a significant volume of fluids injected into the domain.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy, Jutul</span></span>
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
<span class="line"><span>│ Linear solver  │    10.252 │       8.94326 │ 1261 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   20.5041 │       17.8865 │ 2522 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2382 │     5.40 % │ 0.0974 │</span></span>
<span class="line"><span>│ Equations     │ 0.2101 │     6.40 % │ 0.1155 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2543 │     7.74 % │ 0.1399 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2915 │     6.60 % │ 0.1192 │</span></span>
<span class="line"><span>│ Linear setup  │ 1.8425 │    41.73 % │ 0.7536 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1888 │    26.36 % │ 0.4761 │</span></span>
<span class="line"><span>│ Update        │ 0.0752 │     1.70 % │ 0.0307 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0694 │     2.11 % │ 0.0382 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0391 │     0.31 % │ 0.0055 │</span></span>
<span class="line"><span>│ Other         │ 0.0728 │     1.65 % │ 0.0298 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.4157 │   100.00 % │ 1.8060 │</span></span>
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
<span class="line"><span>│                │ 123 steps │ 166 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   4.86179 │       3.60241 │  598 (0) │</span></span>
<span class="line"><span>│ Linearization  │   6.21138 │       4.60241 │  764 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   22.7317 │       16.8434 │ 2796 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   45.4634 │       33.6867 │ 5592 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2424 │     4.48 % │ 0.1450 │</span></span>
<span class="line"><span>│ Equations     │ 0.2350 │     5.54 % │ 0.1795 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2588 │     6.10 % │ 0.1977 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3951 │     7.29 % │ 0.2363 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2187 │    40.96 % │ 1.3268 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1804 │    31.15 % │ 1.0088 │</span></span>
<span class="line"><span>│ Update        │ 0.0788 │     1.46 % │ 0.0471 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0717 │     1.69 % │ 0.0547 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0355 │     0.18 % │ 0.0059 │</span></span>
<span class="line"><span>│ Other         │ 0.0620 │     1.15 % │ 0.0371 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.4162 │   100.00 % │ 3.2389 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬────────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │      Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 179 ministeps │   (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼────────────┤</span></span>
<span class="line"><span>│ Newton         │   6.06504 │        4.1676 │   746 (75) │</span></span>
<span class="line"><span>│ Linearization  │   7.52033 │        5.1676 │   925 (80) │</span></span>
<span class="line"><span>│ Linear solver  │   13.5854 │        9.3352 │  1671 (75) │</span></span>
<span class="line"><span>│ Precond apply  │   27.1707 │       18.6704 │ 3342 (150) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴────────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2399 │     4.83 % │ 0.1790 │</span></span>
<span class="line"><span>│ Equations     │ 0.2348 │     5.87 % │ 0.2172 │</span></span>
<span class="line"><span>│ Assembly      │ 0.3915 │     9.78 % │ 0.3621 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2437 │     4.91 % │ 0.1818 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.0560 │    41.43 % │ 1.5337 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1841 │    16.62 % │ 0.6152 │</span></span>
<span class="line"><span>│ Update        │ 0.0810 │     1.63 % │ 0.0604 │</span></span>
<span class="line"><span>│ Convergence   │ 0.4978 │    12.44 % │ 0.4604 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0335 │     0.16 % │ 0.0060 │</span></span>
<span class="line"><span>│ Other         │ 0.1157 │     2.33 % │ 0.0863 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.9626 │   100.00 % │ 3.7021 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 148 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   3.61789 │       3.00676 │  445 (0) │</span></span>
<span class="line"><span>│ Linearization  │   4.82114 │       4.00676 │  593 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   16.3333 │       13.5743 │ 2009 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   32.6667 │       27.1486 │ 4018 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2369 │     4.50 % │ 0.1054 │</span></span>
<span class="line"><span>│ Equations     │ 0.2217 │     5.60 % │ 0.1314 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2520 │     6.37 % │ 0.1495 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3678 │     6.98 % │ 0.1637 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.1929 │    41.61 % │ 0.9759 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1775 │    30.41 % │ 0.7132 │</span></span>
<span class="line"><span>│ Update        │ 0.0745 │     1.41 % │ 0.0332 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0664 │     1.68 % │ 0.0394 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0332 │     0.21 % │ 0.0049 │</span></span>
<span class="line"><span>│ Other         │ 0.0646 │     1.23 % │ 0.0287 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.2701 │   100.00 % │ 2.3452 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬─────────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │       Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 351 ministeps │    (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼─────────────┤</span></span>
<span class="line"><span>│ Newton         │   19.3333 │       6.77493 │ 2378 (1245) │</span></span>
<span class="line"><span>│ Linearization  │    22.187 │       7.77493 │ 2729 (1328) │</span></span>
<span class="line"><span>│ Linear solver  │   29.0976 │       10.1966 │ 3579 (1286) │</span></span>
<span class="line"><span>│ Precond apply  │   58.1951 │       20.3932 │ 7158 (2572) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴─────────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2533 │     7.10 % │ 0.6024 │</span></span>
<span class="line"><span>│ Equations     │ 0.2194 │     7.05 % │ 0.5987 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2534 │     8.15 % │ 0.6917 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2089 │     5.85 % │ 0.4968 │</span></span>
<span class="line"><span>│ Linear setup  │ 1.7889 │    50.12 % │ 4.2539 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1810 │    15.26 % │ 1.2954 │</span></span>
<span class="line"><span>│ Update        │ 0.0743 │     2.08 % │ 0.1767 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0696 │     2.24 % │ 0.1899 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0254 │     0.10 % │ 0.0089 │</span></span>
<span class="line"><span>│ Other         │ 0.0724 │     2.03 % │ 0.1722 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 3.5688 │   100.00 % │ 8.4866 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬────────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │      Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 189 ministeps │   (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼────────────┤</span></span>
<span class="line"><span>│ Newton         │   6.70732 │       4.36508 │  825 (105) │</span></span>
<span class="line"><span>│ Linearization  │    8.2439 │       5.36508 │ 1014 (112) │</span></span>
<span class="line"><span>│ Linear solver  │   14.3984 │       9.37037 │ 1771 (105) │</span></span>
<span class="line"><span>│ Precond apply  │   28.7967 │       18.7407 │ 3542 (210) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴────────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2395 │     5.72 % │ 0.1976 │</span></span>
<span class="line"><span>│ Equations     │ 0.2698 │     7.92 % │ 0.2736 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2555 │     7.50 % │ 0.2591 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2413 │     5.76 % │ 0.1991 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.0482 │    48.92 % │ 1.6898 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1847 │    18.94 % │ 0.6540 │</span></span>
<span class="line"><span>│ Update        │ 0.0755 │     1.80 % │ 0.0623 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0696 │     2.04 % │ 0.0706 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0310 │     0.17 % │ 0.0059 │</span></span>
<span class="line"><span>│ Other         │ 0.0509 │     1.22 % │ 0.0420 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.1864 │   100.00 % │ 3.4538 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬───────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │     Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 176 ministeps │  (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼───────────┤</span></span>
<span class="line"><span>│ Newton         │   5.68293 │       3.97159 │  699 (45) │</span></span>
<span class="line"><span>│ Linearization  │   7.11382 │       4.97159 │  875 (48) │</span></span>
<span class="line"><span>│ Linear solver  │   13.5691 │       9.48295 │ 1669 (45) │</span></span>
<span class="line"><span>│ Precond apply  │   27.1382 │       18.9659 │ 3338 (90) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴───────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2368 │     5.47 % │ 0.1655 │</span></span>
<span class="line"><span>│ Equations     │ 0.2242 │     6.49 % │ 0.1962 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2534 │     7.33 % │ 0.2217 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.2541 │     5.87 % │ 0.1776 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.0802 │    48.07 % │ 1.4541 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1829 │    20.19 % │ 0.6107 │</span></span>
<span class="line"><span>│ Update        │ 0.0752 │     1.74 % │ 0.0526 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0692 │     2.00 % │ 0.0605 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.2738 │     1.59 % │ 0.0482 │</span></span>
<span class="line"><span>│ Other         │ 0.0543 │     1.25 % │ 0.0379 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 4.3277 │   100.00 % │ 3.0251 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 152 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   3.93496 │       3.18421 │  484 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.17073 │       4.18421 │  636 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   17.7317 │       14.3487 │ 2181 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   35.4634 │       28.6974 │ 4362 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2374 │     4.51 % │ 0.1149 │</span></span>
<span class="line"><span>│ Equations     │ 0.2227 │     5.56 % │ 0.1417 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2542 │     6.34 % │ 0.1617 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3723 │     7.07 % │ 0.1802 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2029 │    41.82 % │ 1.0662 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1764 │    30.18 % │ 0.7694 │</span></span>
<span class="line"><span>│ Update        │ 0.0754 │     1.43 % │ 0.0365 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0678 │     1.69 % │ 0.0431 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0350 │     0.21 % │ 0.0053 │</span></span>
<span class="line"><span>│ Other         │ 0.0634 │     1.20 % │ 0.0307 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.2677 │   100.00 % │ 2.5496 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 153 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   4.05691 │       3.26144 │  499 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.30081 │       4.26144 │  652 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   16.8293 │       13.5294 │ 2070 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   33.6585 │       27.0588 │ 4140 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2411 │     4.62 % │ 0.1203 │</span></span>
<span class="line"><span>│ Equations     │ 0.2778 │     6.95 % │ 0.1811 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2546 │     6.37 % │ 0.1660 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3502 │     6.71 % │ 0.1748 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2161 │    42.43 % │ 1.1058 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1786 │    28.37 % │ 0.7393 │</span></span>
<span class="line"><span>│ Update        │ 0.0757 │     1.45 % │ 0.0378 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0686 │     1.72 % │ 0.0447 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0341 │     0.20 % │ 0.0052 │</span></span>
<span class="line"><span>│ Other         │ 0.0623 │     1.19 % │ 0.0311 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.2227 │   100.00 % │ 2.6061 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 145 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   3.57724 │       3.03448 │  440 (0) │</span></span>
<span class="line"><span>│ Linearization  │    4.7561 │       4.03448 │  585 (0) │</span></span>
<span class="line"><span>│ Linear solver  │   14.3659 │       12.1862 │ 1767 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   28.7317 │       24.3724 │ 3534 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2386 │     4.58 % │ 0.1050 │</span></span>
<span class="line"><span>│ Equations     │ 0.2272 │     5.80 % │ 0.1329 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2549 │     6.51 % │ 0.1491 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3490 │     6.70 % │ 0.1535 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2122 │    42.50 % │ 0.9734 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1792 │    27.65 % │ 0.6333 │</span></span>
<span class="line"><span>│ Update        │ 0.0774 │     1.49 % │ 0.0341 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0703 │     1.79 % │ 0.0411 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.2590 │     1.64 % │ 0.0376 │</span></span>
<span class="line"><span>│ Other         │ 0.0694 │     1.33 % │ 0.0305 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.2056 │   100.00 % │ 2.2905 │</span></span>
<span class="line"><span>╰───────────────┴────────┴────────────┴────────╯</span></span>
<span class="line"><span>Jutul: Simulating 9 years, 44.69 weeks as 123 report steps</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 123 steps │ 150 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │   4.01626 │       3.29333 │  494 (0) │</span></span>
<span class="line"><span>│ Linearization  │   5.23577 │       4.29333 │  644 (0) │</span></span>
<span class="line"><span>│ Linear solver  │    19.065 │       15.6333 │ 2345 (0) │</span></span>
<span class="line"><span>│ Precond apply  │   38.1301 │       31.2667 │ 4690 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬────────┬────────────┬────────╮</span></span>
<span class="line"><span>│ Timing type   │   Each │   Relative │  Total │</span></span>
<span class="line"><span>│               │     ms │ Percentage │      s │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Properties    │ 0.2401 │     4.39 % │ 0.1186 │</span></span>
<span class="line"><span>│ Equations     │ 0.2296 │     5.47 % │ 0.1478 │</span></span>
<span class="line"><span>│ Assembly      │ 0.2580 │     6.15 % │ 0.1662 │</span></span>
<span class="line"><span>│ Linear solve  │ 0.3983 │     7.28 % │ 0.1967 │</span></span>
<span class="line"><span>│ Linear setup  │ 2.2478 │    41.07 % │ 1.1104 │</span></span>
<span class="line"><span>│ Precond apply │ 0.1788 │    31.01 % │ 0.8384 │</span></span>
<span class="line"><span>│ Update        │ 0.0791 │     1.45 % │ 0.0391 │</span></span>
<span class="line"><span>│ Convergence   │ 0.0728 │     1.73 % │ 0.0469 │</span></span>
<span class="line"><span>│ Input/Output  │ 0.0404 │     0.22 % │ 0.0061 │</span></span>
<span class="line"><span>│ Other         │ 0.0680 │     1.24 % │ 0.0336 │</span></span>
<span class="line"><span>├───────────────┼────────┼────────────┼────────┤</span></span>
<span class="line"><span>│ Total         │ 5.4734 │   100.00 % │ 2.7039 │</span></span>
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
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+r+'" alt=""></p><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/workflow/five_spot_ensemble.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/workflow/five_spot_ensemble.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>This example took 44.468571814 seconds to complete.</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',20))])}const A=e(E,[["render",o]]);export{v as __pageData,A as default};
