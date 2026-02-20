import{_ as a,c as n,a5 as i,o as p}from"./chunks/framework.B4o_rkVA.js";const l="/JutulDarcy.jl/v0.2.33/assets/qjexpwy.C6e-zcDU.jpeg",e="/JutulDarcy.jl/v0.2.33/assets/dpvzqef.DxqjO6oT.jpeg",t="/JutulDarcy.jl/v0.2.33/assets/tdzpstr.eGrWMuaZ.jpeg",h="/JutulDarcy.jl/v0.2.33/assets/tnrqnhy.DEQxvZqt.jpeg",k="/JutulDarcy.jl/v0.2.33/assets/qcwpfve.DlpZtOvj.jpeg",r="/JutulDarcy.jl/v0.2.33/assets/hvxszyg.CDQx-kcn.jpeg",E="/JutulDarcy.jl/v0.2.33/assets/kmnliyw.f1yxBKLS.jpeg",c="/JutulDarcy.jl/v0.2.33/assets/gouexnm.DwA6vf8i.jpeg",d="/JutulDarcy.jl/v0.2.33/assets/tgivlii.R2ZEwr7V.jpeg",g="/JutulDarcy.jl/v0.2.33/assets/zjrgerb.CZ7-Cf4r.jpeg",y="/JutulDarcy.jl/v0.2.33/assets/nllkjcz.DwVewtTV.jpeg",D=JSON.parse('{"title":"Intro to sensitivities in JutulDarcy","description":"","frontmatter":{},"headers":[],"relativePath":"examples/intro_sensitivities.md","filePath":"examples/intro_sensitivities.md","lastUpdated":null}'),o={name:"examples/intro_sensitivities.md"};function u(F,s,b,C,A,m){return p(),n("div",null,s[0]||(s[0]=[i(`<h1 id="Intro-to-sensitivities-in-JutulDarcy" tabindex="-1">Intro to sensitivities in JutulDarcy <a class="header-anchor" href="#Intro-to-sensitivities-in-JutulDarcy" aria-label="Permalink to &quot;Intro to sensitivities in JutulDarcy {#Intro-to-sensitivities-in-JutulDarcy}&quot;">​</a></h1><p>Sensitivites with respect to custom parameters: We demonstrate how to set up a simple conceptual model, add new parameters and variable definitions in the form of a new relative permeability function, and calculate and visualize parameter sensitivities.</p><p>We first set up a quarter-five-spot model where the domain is flooded from left to right. Some cells have lower permeability to impede flow and make the scenario more interesting.</p><p>For more details, see the paper <a href="https://doi.org/10.3997/2214-4609.202437111" target="_blank" rel="noreferrer">JutulDarcy.jl - a Fully Differentiable High-Performance Reservoir Simulator Based on Automatic Differentiation</a>.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Jutul, JutulDarcy, GLMakie, HYPRE</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">darcy, kg, meter, year, day, bar </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> si_units</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:darcy</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:kilogram</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:meter</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:year</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:day</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:bar</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1000.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">meter</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">H </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 100.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">meter</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">big </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> false</span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"> # Paper uses big, takes some more time to run</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> big</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    nx </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 500</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">else</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    nx </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 100</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dx </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> L</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">nx</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">g </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> CartesianMesh</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">((nx, nx, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), (L, L, H))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">nc </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> number_of_cells</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(g)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">perm </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> fill</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">darcy, nc)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">reservoir </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reservoir_domain</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(g, permeability </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">darcy)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">centroids </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> reservoir[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:cell_centroids</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rock_type </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> fill</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, nc)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (i, x, y) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> zip</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">eachindex</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(perm), centroids[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :], centroids[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    xseg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.75</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    yseg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.75</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> xseg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">||</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> yseg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        rock_type[i] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 2</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    xseg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.55</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.50</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.55</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    yseg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.55</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.50</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.55</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> xseg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">||</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> yseg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        rock_type[i] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 3</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    xseg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.25</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    yseg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.25</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">L)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> xseg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">||</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> yseg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        rock_type[i] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 4</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">perm </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> reservoir[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:permeability</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> perm[rock_type </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">darcy</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> perm[rock_type </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.005</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">darcy</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> perm[rock_type </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.01</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">darcy</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">I </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_vertical_well</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(reservoir, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, name </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Injector</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">P </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_vertical_well</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(reservoir, nx, nx, name </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Producer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">phases </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AqueousPhase</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">VaporPhase</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rhoWS, rhoGS </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1000.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">kg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">meter</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">700.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">kg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">meter</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">system </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ImmiscibleSystem</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(phases, reference_densities </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (rhoWS, rhoGS))</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">model, </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_reservoir_model</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(reservoir, system, wells </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [I, P])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rmodel </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reservoir_model</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>SimulationModel:</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Model with 20000 degrees of freedom, 20000 equations and 79600 parameters</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  domain:</span></span>
<span class="line"><span>    DiscretizedDomain with MinimalTPFATopology (10000 cells, 19800 faces) and discretizations for mass_flow, heat_flow</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  system:</span></span>
<span class="line"><span>    ImmiscibleSystem with AqueousPhase, VaporPhase</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  context:</span></span>
<span class="line"><span>    DefaultContext(BlockMajorLayout(false), 1000, 1)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  formulation:</span></span>
<span class="line"><span>    FullyImplicitFormulation()</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  data_domain:</span></span>
<span class="line"><span>    DataDomain wrapping CartesianMesh (3D) with 100x100x1=10000 cells with 17 data fields added:</span></span>
<span class="line"><span>  10000 Cells</span></span>
<span class="line"><span>    :permeability =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :porosity =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :rock_thermal_conductivity =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :fluid_thermal_conductivity =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :rock_density =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :cell_centroids =&gt; 3×10000 Matrix{Float64}</span></span>
<span class="line"><span>    :volumes =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>  19800 Faces</span></span>
<span class="line"><span>    :neighbors =&gt; 2×19800 Matrix{Int64}</span></span>
<span class="line"><span>    :areas =&gt; 19800 Vector{Float64}</span></span>
<span class="line"><span>    :normals =&gt; 3×19800 Matrix{Float64}</span></span>
<span class="line"><span>    :face_centroids =&gt; 3×19800 Matrix{Float64}</span></span>
<span class="line"><span>  39600 HalfFaces</span></span>
<span class="line"><span>    :half_face_cells =&gt; 39600 Vector{Int64}</span></span>
<span class="line"><span>    :half_face_faces =&gt; 39600 Vector{Int64}</span></span>
<span class="line"><span>  20400 BoundaryFaces</span></span>
<span class="line"><span>    :boundary_areas =&gt; 20400 Vector{Float64}</span></span>
<span class="line"><span>    :boundary_centroids =&gt; 3×20400 Matrix{Float64}</span></span>
<span class="line"><span>    :boundary_normals =&gt; 3×20400 Matrix{Float64}</span></span>
<span class="line"><span>    :boundary_neighbors =&gt; 20400 Vector{Int64}</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  primary_variables:</span></span>
<span class="line"><span>   1) Pressure    ∈ 10000 Cells: 1 dof each</span></span>
<span class="line"><span>   2) Saturations ∈ 10000 Cells: 1 dof, 2 values each</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  secondary_variables:</span></span>
<span class="line"><span>   1) PhaseMassDensities     ∈ 10000 Cells: 2 values each</span></span>
<span class="line"><span>      -&gt; ConstantCompressibilityDensities as evaluator</span></span>
<span class="line"><span>   2) TotalMasses            ∈ 10000 Cells: 2 values each</span></span>
<span class="line"><span>      -&gt; TotalMasses as evaluator</span></span>
<span class="line"><span>   3) RelativePermeabilities ∈ 10000 Cells: 2 values each</span></span>
<span class="line"><span>      -&gt; BrooksCoreyRelativePermeabilities as evaluator</span></span>
<span class="line"><span>   4) PhaseMobilities        ∈ 10000 Cells: 2 values each</span></span>
<span class="line"><span>      -&gt; JutulDarcy.PhaseMobilities as evaluator</span></span>
<span class="line"><span>   5) PhaseMassMobilities    ∈ 10000 Cells: 2 values each</span></span>
<span class="line"><span>      -&gt; JutulDarcy.PhaseMassMobilities as evaluator</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  parameters:</span></span>
<span class="line"><span>   1) Transmissibilities        ∈ 19800 Faces: Scalar</span></span>
<span class="line"><span>   2) TwoPointGravityDifference ∈ 19800 Faces: Scalar</span></span>
<span class="line"><span>   3) ConnateWater              ∈ 10000 Cells: Scalar</span></span>
<span class="line"><span>   4) PhaseViscosities          ∈ 10000 Cells: 2 values each</span></span>
<span class="line"><span>   5) FluidVolume               ∈ 10000 Cells: Scalar</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  equations:</span></span>
<span class="line"><span>   1) mass_conservation ∈ 10000 Cells: 2 values each</span></span>
<span class="line"><span>      -&gt; ConservationLaw{:TotalMasses, TwoPointPotentialFlowHardCoded{Vector{Int64}, Vector{@NamedTuple{self::Int64, other::Int64, face::Int64, face_sign::Int64}}}, Jutul.DefaultFlux, 2}</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  output_variables:</span></span>
<span class="line"><span>    Pressure, Saturations, TotalMasses</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  extra:</span></span>
<span class="line"><span>    OrderedCollections.OrderedDict{Symbol, Any} with keys: Symbol[]</span></span></code></pre></div><h2 id="Plot-the-initial-variable-graph" tabindex="-1">Plot the initial variable graph <a class="header-anchor" href="#Plot-the-initial-variable-graph" aria-label="Permalink to &quot;Plot the initial variable graph {#Plot-the-initial-variable-graph}&quot;">​</a></h2><p>We plot the default variable graph that describes how the different variables relate to each other. When we add a new parameter and property in the next section, the graph is automatically modified.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> NetworkLayout, LayeredLayouts, GraphMakie</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Jutul</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_variable_graph</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rmodel)</span></span></code></pre></div><p><img src="`+l+`" alt=""></p><h2 id="Change-the-variables" tabindex="-1">Change the variables <a class="header-anchor" href="#Change-the-variables" aria-label="Permalink to &quot;Change the variables {#Change-the-variables}&quot;">​</a></h2><p>We replace the density variable with a more compressible version, and we also define a new relative permeability variable that depends on a new parameter <code>KrExponents</code> to define the exponent of the relative permeability in each cell and phase of the model.</p><p>This is done through several steps:</p><ol><li><p>First, we define the type</p></li><li><p>We define functions that act on that type, in particular the update function that is used to evaluate the new relative permeability during the simulation for named inputs <code>Saturations</code> and <code>KrExponents</code>.</p></li><li><p>We define the <code>KrExponents</code> as a model parameter with a default value, that can subsequently be used by the relative permeability.</p></li></ol><p>Finally we plot the variable graph again to verify that the new relationship has been included in our model.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">c </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-6</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-4</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">density </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ConstantCompressibilityDensities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(p_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar, density_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [rhoWS, rhoGS], compressibility </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> c)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">replace_variables!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rmodel, PhaseMassDensities </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> density);</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">import</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> AbstractRelativePermeabilities, PhaseVariables</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> MyKr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> AbstractRelativePermeabilities</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> end</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@jutul_secondary</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> update_my_kr!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(vals, def</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">MyKr</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, model, Saturations, KrExponents, cells_to_update)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> c </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cells_to_update</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ph </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> axes</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(vals, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            S_α </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> max</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(Saturations[ph, c], </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            n_α </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> KrExponents[ph, c]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            vals[ph, c] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> S_α</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">n_α</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> MyKrExp </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> PhaseVariables</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Jutul</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">default_value</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">MyKrExp</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 2.0</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">set_parameters!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rmodel, KrExponents </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> MyKrExp</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">replace_variables!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rmodel, RelativePermeabilities </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> MyKr</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">());</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Jutul</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_variable_graph</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rmodel)</span></span></code></pre></div><p><img src="`+e+`" alt=""></p><h2 id="Set-up-scenario-and-simulate" tabindex="-1">Set up scenario and simulate <a class="header-anchor" href="#Set-up-scenario-and-simulate" aria-label="Permalink to &quot;Set up scenario and simulate {#Set-up-scenario-and-simulate}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">parameters </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_parameters</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">exponents </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> parameters[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:KrExponents</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (cell, rtype) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> enumerate</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rock_type)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> rtype </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        exp_w </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 2</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        exp_g </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 3</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    else</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        exp_w </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        exp_g </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 2</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    exponents[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, cell] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> exp_w</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    exponents[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, cell] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> exp_g</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">pv </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> pore_volume</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, parameters)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">state0 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_reservoir_state</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, Pressure </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 150</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar, Saturations </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> repeat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">30.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">day, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">12</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">pv </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> pore_volume</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, parameters)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">total_time </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sum</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(dt)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">inj_rate </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sum</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(pv)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">total_time</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rate_target </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> TotalRateTarget</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(inj_rate)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">I_ctrl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> InjectorControl</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rate_target, [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], density </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> rhoGS)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bhp_target </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> BottomHolePressureTarget</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">50</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bar)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">P_ctrl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ProducerControl</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(bhp_target)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">controls </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Dict</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">controls[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Injector</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> I_ctrl</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">controls[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Producer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> P_ctrl</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">forces </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_reservoir_forces</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, control </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> controls)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">case </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> JutulCase</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, dt, forces, parameters </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> parameters, state0 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> state0)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">result </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case, output_substates </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">);</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ws, states </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> result</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ws</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Producer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:grat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span></span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps   3%| |  ETA: 0:04:15\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 2/60 (3.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     56 iterations in 7.82 s (139.60 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps   5%| |  ETA: 0:02:55\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 3/60 (5.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     73 iterations in 8.11 s (111.16 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps   7%|▏|  ETA: 0:02:12\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 4/60 (6.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     85 iterations in 8.33 s (97.97 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps   8%|▏|  ETA: 0:01:46\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 5/60 (8.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     96 iterations in 8.52 s (88.79 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  10%|▏|  ETA: 0:01:28\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 6/60 (10.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     106 iterations in 8.70 s (82.04 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  11%|▏|  ETA: 0:01:16\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 7/60 (11.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     116 iterations in 8.87 s (76.47 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  13%|▏|  ETA: 0:01:06\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 8/60 (13.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     125 iterations in 9.02 s (72.17 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  15%|▏|  ETA: 0:00:58\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 9/60 (15.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     135 iterations in 9.19 s (68.06 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  16%|▏|  ETA: 0:00:52\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 10/60 (16.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     144 iterations in 9.33 s (64.80 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  18%|▏|  ETA: 0:00:47\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 11/60 (18.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     152 iterations in 9.47 s (62.28 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  20%|▎|  ETA: 0:00:43\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 12/60 (20.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     161 iterations in 9.64 s (59.88 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  21%|▎|  ETA: 0:00:40\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 13/60 (21.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     170 iterations in 9.80 s (57.64 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  23%|▎|  ETA: 0:00:37\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 14/60 (23.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     180 iterations in 9.97 s (55.40 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  25%|▎|  ETA: 0:00:34\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 15/60 (25.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     189 iterations in 10.11 s (53.50 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  26%|▎|  ETA: 0:00:31\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 16/60 (26.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     197 iterations in 10.24 s (51.97 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  28%|▎|  ETA: 0:00:29\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 17/60 (28.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     205 iterations in 10.36 s (50.56 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  30%|▎|  ETA: 0:00:27\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 18/60 (30.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     213 iterations in 10.49 s (49.24 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  31%|▎|  ETA: 0:00:26\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 19/60 (31.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     222 iterations in 10.66 s (48.01 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  33%|▍|  ETA: 0:00:24\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 20/60 (33.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     230 iterations in 10.79 s (46.91 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  34%|▍|  ETA: 0:00:23\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 21/60 (35.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     239 iterations in 10.94 s (45.77 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  36%|▍|  ETA: 0:00:21\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 22/60 (36.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     249 iterations in 11.11 s (44.62 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  38%|▍|  ETA: 0:00:20\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 23/60 (38.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     258 iterations in 11.27 s (43.68 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  39%|▍|  ETA: 0:00:19\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 24/60 (40.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     266 iterations in 11.40 s (42.87 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  41%|▍|  ETA: 0:00:18\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 25/60 (41.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     274 iterations in 11.53 s (42.09 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  43%|▍|  ETA: 0:00:17\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 26/60 (43.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     283 iterations in 11.72 s (41.41 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  44%|▌|  ETA: 0:00:16\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 27/60 (45.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     291 iterations in 11.86 s (40.74 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  46%|▌|  ETA: 0:00:15\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 28/60 (46.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     299 iterations in 11.99 s (40.09 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  48%|▌|  ETA: 0:00:14\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 29/60 (48.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     307 iterations in 12.13 s (39.52 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  49%|▌|  ETA: 0:00:14\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 30/60 (50.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     315 iterations in 12.28 s (38.99 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  51%|▌|  ETA: 0:00:13\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 31/60 (51.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     323 iterations in 12.41 s (38.43 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  52%|▌|  ETA: 0:00:12\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 32/60 (53.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     331 iterations in 12.55 s (37.91 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  54%|▌|  ETA: 0:00:12\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 33/60 (55.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     339 iterations in 12.69 s (37.42 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  56%|▌|  ETA: 0:00:11\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 34/60 (56.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     347 iterations in 12.85 s (37.03 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  57%|▋|  ETA: 0:00:10\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 35/60 (58.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     355 iterations in 12.99 s (36.59 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  59%|▋|  ETA: 0:00:10\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 36/60 (60.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     364 iterations in 13.14 s (36.10 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  61%|▋|  ETA: 0:00:09\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 37/60 (61.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     372 iterations in 13.28 s (35.70 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  62%|▋|  ETA: 0:00:09\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 38/60 (63.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     381 iterations in 13.45 s (35.29 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  64%|▋|  ETA: 0:00:08\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 39/60 (65.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     390 iterations in 13.61 s (34.89 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  66%|▋|  ETA: 0:00:08\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 40/60 (66.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     399 iterations in 13.76 s (34.50 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  67%|▋|  ETA: 0:00:07\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 41/60 (68.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     408 iterations in 13.96 s (34.22 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  69%|▊|  ETA: 0:00:07\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 42/60 (70.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     418 iterations in 14.13 s (33.81 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  70%|▊|  ETA: 0:00:06\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 43/60 (71.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     428 iterations in 14.30 s (33.41 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  72%|▊|  ETA: 0:00:06\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 44/60 (73.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     438 iterations in 14.47 s (33.04 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  74%|▊|  ETA: 0:00:06\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 45/60 (75.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     449 iterations in 14.66 s (32.65 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  75%|▊|  ETA: 0:00:05\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 46/60 (76.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     466 iterations in 14.98 s (32.16 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  77%|▊|  ETA: 0:00:05\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 47/60 (78.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     477 iterations in 15.20 s (31.86 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  79%|▊|  ETA: 0:00:04\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 48/60 (80.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     483 iterations in 15.31 s (31.69 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  82%|▉|  ETA: 0:00:04\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 50/60 (83.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     491 iterations in 15.45 s (31.47 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  85%|▉|  ETA: 0:00:03\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 52/60 (86.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     499 iterations in 15.60 s (31.27 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  89%|▉|  ETA: 0:00:02\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 54/60 (90.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     507 iterations in 15.74 s (31.05 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  92%|▉|  ETA: 0:00:02\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 56/60 (93.33% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     515 iterations in 15.89 s (30.85 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  95%|█|  ETA: 0:00:01\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 58/60 (96.67% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     523 iterations in 16.05 s (30.69 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps  98%|█|  ETA: 0:00:00\x1B[K</span></span>
<span class="line"><span>  Progress:  Solving step 60/60 (100.00% of time interval complete)\x1B[K</span></span>
<span class="line"><span>  Stats:     531 iterations in 16.19 s (30.50 ms each)\x1B[K</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span>\x1B[A</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>\x1B[K\x1B[A</span></span>
<span class="line"><span>Simulating 4 years, 48.43 weeks as 60 report steps 100%|█| Time: 0:00:17\x1B[K</span></span>
<span class="line"><span>  Progress:  Solved step 60/60\x1B[K</span></span>
<span class="line"><span>  Stats:     535 iterations in 16.26 s (30.40 ms each)\x1B[K</span></span>
<span class="line"><span>╭────────────────┬──────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │ Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 60 steps │ 117 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼──────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │  8.91667 │       4.57265 │  535 (0) │</span></span>
<span class="line"><span>│ Linearization  │  10.8667 │       5.57265 │  652 (0) │</span></span>
<span class="line"><span>│ Linear solver  │    30.85 │       15.8205 │ 1851 (0) │</span></span>
<span class="line"><span>│ Precond apply  │     61.7 │        31.641 │ 3702 (0) │</span></span>
<span class="line"><span>╰────────────────┴──────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬─────────┬────────────┬─────────╮</span></span>
<span class="line"><span>│ Timing type   │    Each │   Relative │   Total │</span></span>
<span class="line"><span>│               │      ms │ Percentage │       s │</span></span>
<span class="line"><span>├───────────────┼─────────┼────────────┼─────────┤</span></span>
<span class="line"><span>│ Properties    │  0.5998 │     1.97 % │  0.3209 │</span></span>
<span class="line"><span>│ Equations     │  5.3372 │    21.40 % │  3.4799 │</span></span>
<span class="line"><span>│ Assembly      │  3.1869 │    12.78 % │  2.0778 │</span></span>
<span class="line"><span>│ Linear solve  │  1.1479 │     3.78 % │  0.6141 │</span></span>
<span class="line"><span>│ Linear setup  │  8.2018 │    26.98 % │  4.3879 │</span></span>
<span class="line"><span>│ Precond apply │  0.8325 │    18.95 % │  3.0819 │</span></span>
<span class="line"><span>│ Update        │  0.6551 │     2.15 % │  0.3505 │</span></span>
<span class="line"><span>│ Convergence   │  1.3209 │     5.30 % │  0.8612 │</span></span>
<span class="line"><span>│ Input/Output  │  0.4132 │     0.30 % │  0.0483 │</span></span>
<span class="line"><span>│ Other         │  1.9449 │     6.40 % │  1.0405 │</span></span>
<span class="line"><span>├───────────────┼─────────┼────────────┼─────────┤</span></span>
<span class="line"><span>│ Total         │ 30.3983 │   100.00 % │ 16.2631 │</span></span>
<span class="line"><span>╰───────────────┴─────────┴────────────┴─────────╯</span></span>
<span class="line"><span>Legend</span></span>
<span class="line"><span>┌───────┬──────────────────┬──────┬─────────────────────────┐</span></span>
<span class="line"><span>│ Label │ Description      │ Unit │ Type of quantity        │</span></span>
<span class="line"><span>├───────┼──────────────────┼──────┼─────────────────────────┤</span></span>
<span class="line"><span>│ grat  │ Surface gas rate │ m³/s │ surface_volume_per_time │</span></span>
<span class="line"><span>└───────┴──────────────────┴──────┴─────────────────────────┘</span></span>
<span class="line"><span>Producer result</span></span>
<span class="line"><span>┌─────────┬──────────────┐</span></span>
<span class="line"><span>│    time │         grat │</span></span>
<span class="line"><span>│    days │         m³/s │</span></span>
<span class="line"><span>├─────────┼──────────────┤</span></span>
<span class="line"><span>│     1.0 │         -0.0 │</span></span>
<span class="line"><span>│ 1.64286 │         -0.0 │</span></span>
<span class="line"><span>│ 2.05612 │         -0.0 │</span></span>
<span class="line"><span>│     2.8 │         -0.0 │</span></span>
<span class="line"><span>│ 4.47372 │         -0.0 │</span></span>
<span class="line"><span>│ 7.48643 │         -0.0 │</span></span>
<span class="line"><span>│ 12.0055 │         -0.0 │</span></span>
<span class="line"><span>│ 18.7841 │         -0.0 │</span></span>
<span class="line"><span>│  24.392 │         -0.0 │</span></span>
<span class="line"><span>│    30.0 │         -0.0 │</span></span>
<span class="line"><span>│ 38.4119 │         -0.0 │</span></span>
<span class="line"><span>│  49.206 │         -0.0 │</span></span>
<span class="line"><span>│    60.0 │         -0.0 │</span></span>
<span class="line"><span>│  73.878 │         -0.0 │</span></span>
<span class="line"><span>│    90.0 │         -0.0 │</span></span>
<span class="line"><span>│   105.0 │         -0.0 │</span></span>
<span class="line"><span>│   120.0 │         -0.0 │</span></span>
<span class="line"><span>│   135.0 │         -0.0 │</span></span>
<span class="line"><span>│   150.0 │         -0.0 │</span></span>
<span class="line"><span>│   165.0 │         -0.0 │</span></span>
<span class="line"><span>│   180.0 │         -0.0 │</span></span>
<span class="line"><span>│   195.0 │         -0.0 │</span></span>
<span class="line"><span>│   210.0 │         -0.0 │</span></span>
<span class="line"><span>│   225.0 │         -0.0 │</span></span>
<span class="line"><span>│   240.0 │         -0.0 │</span></span>
<span class="line"><span>│   255.0 │         -0.0 │</span></span>
<span class="line"><span>│   270.0 │         -0.0 │</span></span>
<span class="line"><span>│   285.0 │         -0.0 │</span></span>
<span class="line"><span>│   300.0 │         -0.0 │</span></span>
<span class="line"><span>│   315.0 │         -0.0 │</span></span>
<span class="line"><span>│   330.0 │         -0.0 │</span></span>
<span class="line"><span>│   345.0 │         -0.0 │</span></span>
<span class="line"><span>│   360.0 │         -0.0 │</span></span>
<span class="line"><span>│   375.0 │         -0.0 │</span></span>
<span class="line"><span>│   390.0 │         -0.0 │</span></span>
<span class="line"><span>│   405.0 │         -0.0 │</span></span>
<span class="line"><span>│   420.0 │         -0.0 │</span></span>
<span class="line"><span>│   435.0 │         -0.0 │</span></span>
<span class="line"><span>│   450.0 │         -0.0 │</span></span>
<span class="line"><span>│   465.0 │         -0.0 │</span></span>
<span class="line"><span>│   480.0 │         -0.0 │</span></span>
<span class="line"><span>│   495.0 │         -0.0 │</span></span>
<span class="line"><span>│   510.0 │         -0.0 │</span></span>
<span class="line"><span>│   525.0 │         -0.0 │</span></span>
<span class="line"><span>│   540.0 │         -0.0 │</span></span>
<span class="line"><span>│   555.0 │         -0.0 │</span></span>
<span class="line"><span>│   570.0 │         -0.0 │</span></span>
<span class="line"><span>│   585.0 │         -0.0 │</span></span>
<span class="line"><span>│   600.0 │         -0.0 │</span></span>
<span class="line"><span>│   615.0 │         -0.0 │</span></span>
<span class="line"><span>│   630.0 │         -0.0 │</span></span>
<span class="line"><span>│   645.0 │         -0.0 │</span></span>
<span class="line"><span>│   660.0 │         -0.0 │</span></span>
<span class="line"><span>│   675.0 │         -0.0 │</span></span>
<span class="line"><span>│   690.0 │         -0.0 │</span></span>
<span class="line"><span>│   705.0 │         -0.0 │</span></span>
<span class="line"><span>│   720.0 │         -0.0 │</span></span>
<span class="line"><span>│   735.0 │         -0.0 │</span></span>
<span class="line"><span>│   750.0 │         -0.0 │</span></span>
<span class="line"><span>│   765.0 │         -0.0 │</span></span>
<span class="line"><span>│   780.0 │         -0.0 │</span></span>
<span class="line"><span>│   795.0 │         -0.0 │</span></span>
<span class="line"><span>│   810.0 │         -0.0 │</span></span>
<span class="line"><span>│   825.0 │         -0.0 │</span></span>
<span class="line"><span>│   840.0 │         -0.0 │</span></span>
<span class="line"><span>│   855.0 │         -0.0 │</span></span>
<span class="line"><span>│   870.0 │         -0.0 │</span></span>
<span class="line"><span>│   885.0 │         -0.0 │</span></span>
<span class="line"><span>│   900.0 │         -0.0 │</span></span>
<span class="line"><span>│   915.0 │         -0.0 │</span></span>
<span class="line"><span>│   930.0 │         -0.0 │</span></span>
<span class="line"><span>│   945.0 │         -0.0 │</span></span>
<span class="line"><span>│   960.0 │         -0.0 │</span></span>
<span class="line"><span>│   975.0 │         -0.0 │</span></span>
<span class="line"><span>│   990.0 │         -0.0 │</span></span>
<span class="line"><span>│  1005.0 │         -0.0 │</span></span>
<span class="line"><span>│  1020.0 │         -0.0 │</span></span>
<span class="line"><span>│  1035.0 │         -0.0 │</span></span>
<span class="line"><span>│  1050.0 │         -0.0 │</span></span>
<span class="line"><span>│  1065.0 │         -0.0 │</span></span>
<span class="line"><span>│  1080.0 │         -0.0 │</span></span>
<span class="line"><span>│  1095.0 │         -0.0 │</span></span>
<span class="line"><span>│  1110.0 │         -0.0 │</span></span>
<span class="line"><span>│  1125.0 │         -0.0 │</span></span>
<span class="line"><span>│  1140.0 │         -0.0 │</span></span>
<span class="line"><span>│  1155.0 │         -0.0 │</span></span>
<span class="line"><span>│  1170.0 │         -0.0 │</span></span>
<span class="line"><span>│  1185.0 │         -0.0 │</span></span>
<span class="line"><span>│  1200.0 │         -0.0 │</span></span>
<span class="line"><span>│  1215.0 │         -0.0 │</span></span>
<span class="line"><span>│  1230.0 │         -0.0 │</span></span>
<span class="line"><span>│  1245.0 │         -0.0 │</span></span>
<span class="line"><span>│  1260.0 │         -0.0 │</span></span>
<span class="line"><span>│  1275.0 │         -0.0 │</span></span>
<span class="line"><span>│  1290.0 │         -0.0 │</span></span>
<span class="line"><span>│  1305.0 │         -0.0 │</span></span>
<span class="line"><span>│  1320.0 │         -0.0 │</span></span>
<span class="line"><span>│  1335.0 │         -0.0 │</span></span>
<span class="line"><span>│  1350.0 │ -0.000264065 │</span></span>
<span class="line"><span>│  1365.0 │  -0.00945487 │</span></span>
<span class="line"><span>│  1372.5 │   -0.0156987 │</span></span>
<span class="line"><span>│  1380.0 │    -0.020066 │</span></span>
<span class="line"><span>│  1395.0 │   -0.0252017 │</span></span>
<span class="line"><span>│  1410.0 │   -0.0286516 │</span></span>
<span class="line"><span>│  1440.0 │   -0.0325127 │</span></span>
<span class="line"><span>│  1470.0 │   -0.0348996 │</span></span>
<span class="line"><span>│  1500.0 │   -0.0365194 │</span></span>
<span class="line"><span>│  1530.0 │   -0.0377322 │</span></span>
<span class="line"><span>│  1560.0 │   -0.0387128 │</span></span>
<span class="line"><span>│  1590.0 │   -0.0395447 │</span></span>
<span class="line"><span>│  1620.0 │    -0.040269 │</span></span>
<span class="line"><span>│  1650.0 │   -0.0409083 │</span></span>
<span class="line"><span>│  1680.0 │   -0.0414781 │</span></span>
<span class="line"><span>│  1710.0 │   -0.0419904 │</span></span>
<span class="line"><span>│  1740.0 │   -0.0424565 │</span></span>
<span class="line"><span>│  1770.0 │   -0.0428874 │</span></span>
<span class="line"><span>│  1800.0 │   -0.0432951 │</span></span>
<span class="line"><span>└─────────┴──────────────┘</span></span></code></pre></div><h2 id="Define-objective-function" tabindex="-1">Define objective function <a class="header-anchor" href="#Define-objective-function" aria-label="Permalink to &quot;Define objective function {#Define-objective-function}&quot;">​</a></h2><p>We let the objective function be the amount produced of produced gas, normalized by the injected amount.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GLMakie</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> objective_function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, state, Δt, step_i, forces)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    grat </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">compute_well_qoi</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, state, forces, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Producer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, SurfaceGasRateTarget)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Δt</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">grat</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(inj_rate</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">total_time)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">data_domain_with_gradients </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">reservoir_sensitivities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case, result, objective_function, include_parameters </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>DataDomain wrapping CartesianMesh (3D) with 100x100x1=10000 cells with 23 data fields added:</span></span>
<span class="line"><span>  10000 Cells</span></span>
<span class="line"><span>    :permeability =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :porosity =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :rock_thermal_conductivity =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :fluid_thermal_conductivity =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :rock_density =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :cell_centroids =&gt; 3×10000 Matrix{Float64}</span></span>
<span class="line"><span>    :volumes =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :ConnateWater =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :PhaseViscosities =&gt; 2×10000 Matrix{Float64}</span></span>
<span class="line"><span>    :FluidVolume =&gt; 10000 Vector{Float64}</span></span>
<span class="line"><span>    :KrExponents =&gt; 2×10000 Matrix{Float64}</span></span>
<span class="line"><span>  19800 Faces</span></span>
<span class="line"><span>    :neighbors =&gt; 2×19800 Matrix{Int64}</span></span>
<span class="line"><span>    :areas =&gt; 19800 Vector{Float64}</span></span>
<span class="line"><span>    :normals =&gt; 3×19800 Matrix{Float64}</span></span>
<span class="line"><span>    :face_centroids =&gt; 3×19800 Matrix{Float64}</span></span>
<span class="line"><span>    :Transmissibilities =&gt; 19800 Vector{Float64}</span></span>
<span class="line"><span>    :TwoPointGravityDifference =&gt; 19800 Vector{Float64}</span></span>
<span class="line"><span>  39600 HalfFaces</span></span>
<span class="line"><span>    :half_face_cells =&gt; 39600 Vector{Int64}</span></span>
<span class="line"><span>    :half_face_faces =&gt; 39600 Vector{Int64}</span></span>
<span class="line"><span>  20400 BoundaryFaces</span></span>
<span class="line"><span>    :boundary_areas =&gt; 20400 Vector{Float64}</span></span>
<span class="line"><span>    :boundary_centroids =&gt; 3×20400 Matrix{Float64}</span></span>
<span class="line"><span>    :boundary_normals =&gt; 3×20400 Matrix{Float64}</span></span>
<span class="line"><span>    :boundary_neighbors =&gt; 20400 Vector{Int64}</span></span></code></pre></div><h2 id="Launch-interactive-plotter-for-cell-wise-gradients" tabindex="-1">Launch interactive plotter for cell-wise gradients <a class="header-anchor" href="#Launch-interactive-plotter-for-cell-wise-gradients" aria-label="Permalink to &quot;Launch interactive plotter for cell-wise gradients {#Launch-interactive-plotter-for-cell-wise-gradients}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(data_domain_with_gradients)</span></span></code></pre></div><p><img src="`+t+`" alt=""></p><h2 id="Set-up-plotting-functions" tabindex="-1">Set up plotting functions <a class="header-anchor" href="#Set-up-plotting-functions" aria-label="Permalink to &quot;Set up plotting functions {#Set-up-plotting-functions}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">∂K </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data_domain_with_gradients[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:permeability</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">∂ϕ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data_domain_with_gradients[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:porosity</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> get_cscale</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    minv0, maxv0 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> extrema</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    minv </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> min</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(minv0, </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">maxv0)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    maxv </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> max</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(maxv0, </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">minv0)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (minv, maxv)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(title, vals; kwarg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    myplot!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, title, vals; kwarg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> fig</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> myplot!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig, I, J, title, vals; is_grad </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, is_log </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> missing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, contourplot </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, nticks </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, ticks </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> missing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colorbar </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, kwarg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[I, J], title </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> title)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> is_grad</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        if</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ismissing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(colorrange)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> get_cscale</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(vals)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        cmap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :seismic</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    else</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        if</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ismissing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(colorrange)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> extrema</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(vals)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        cmap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :seaborn_icefire_gradient</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    hidedecorations!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    hidespines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    arg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (; colormap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cmap, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> colorrange, kwarg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    plt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> plot_cell_data!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, g, vals; shading </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> NoShading, arg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> colorbar</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        if</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ismissing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ticks)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            ticks </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> range</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(colorrange</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, nticks)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        Colorbar</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[I, J</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], plt, ticks </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ticks, ticklabelsize </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 25</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 25</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> fig</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>myplot! (generic function with 1 method)</span></span></code></pre></div><h2 id="Plot-the-permeability" tabindex="-1">Plot the permeability <a class="header-anchor" href="#Plot-the-permeability" aria-label="Permalink to &quot;Plot the permeability {#Plot-the-permeability}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Permeability&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, perm</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">./</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">darcy, colorscale </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> log10, ticks </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.01</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span></code></pre></div><p><img src="`+h+`" alt=""></p><h2 id="Plot-the-evolution-of-the-gas-saturation" tabindex="-1">Plot the evolution of the gas saturation <a class="header-anchor" href="#Plot-the-evolution-of-the-gas-saturation" aria-label="Permalink to &quot;Plot the evolution of the gas saturation {#Plot-the-evolution-of-the-gas-saturation}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1200</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">400</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> states[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">25</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Saturations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Gas saturation&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, sg, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), colorbar </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> states[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">70</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Saturations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Gas saturation&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, sg, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), colorbar </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> states[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">end</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Saturations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Gas saturation&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, sg, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># ## Plot the sensitivity of the objective with respect to permeability</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> big</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cticks </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0005</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0005</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">else</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.05</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.05</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cticks </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.05</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.025</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.025</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.05</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;perm_sens&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, ∂K</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">darcy, is_grad </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, ticks </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cticks, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cr)</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># ## Plot the sensitivity of the objective with respect to porosity</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> big</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.00001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.00001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">else</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.00025</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.00025</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;porosity_sens&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, ∂ϕ, is_grad </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cr)</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">∂xyz </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data_domain_with_gradients[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:cell_centroids</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">∂x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ∂xyz[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">∂y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ∂xyz[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">∂z </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ∂xyz[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> big</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">else</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-7</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-7</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>2-element Vector{Float64}:</span></span>
<span class="line"><span> -1.0e-7</span></span>
<span class="line"><span>  1.0e-7</span></span></code></pre></div><h2 id="Plot-the-sensitivity-of-the-objective-with-respect-to-x-cell-centroids" tabindex="-1">Plot the sensitivity of the objective with respect to x cell centroids <a class="header-anchor" href="#Plot-the-sensitivity-of-the-objective-with-respect-to-x-cell-centroids" aria-label="Permalink to &quot;Plot the sensitivity of the objective with respect to x cell centroids {#Plot-the-sensitivity-of-the-objective-with-respect-to-x-cell-centroids}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;dx_sens&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, ∂x, is_grad </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cr)</span></span></code></pre></div><p><img src="`+k+'" alt=""></p><h2 id="Plot-the-sensitivity-of-the-objective-with-respect-to-y-cell-centroids" tabindex="-1">Plot the sensitivity of the objective with respect to y cell centroids <a class="header-anchor" href="#Plot-the-sensitivity-of-the-objective-with-respect-to-y-cell-centroids" aria-label="Permalink to &quot;Plot the sensitivity of the objective with respect to y cell centroids {#Plot-the-sensitivity-of-the-objective-with-respect-to-y-cell-centroids}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;dy_sens&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, ∂y, is_grad </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cr)</span></span></code></pre></div><p><img src="'+r+'" alt=""></p><h2 id="Plot-the-sensitivity-of-the-objective-with-respect-to-z-cell-centroids" tabindex="-1">Plot the sensitivity of the objective with respect to z cell centroids <a class="header-anchor" href="#Plot-the-sensitivity-of-the-objective-with-respect-to-z-cell-centroids" aria-label="Permalink to &quot;Plot the sensitivity of the objective with respect to z cell centroids {#Plot-the-sensitivity-of-the-objective-with-respect-to-z-cell-centroids}&quot;">​</a></h2><p>Note: The effect here is primarily coming from gravity.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;dz_sens&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, ∂z, is_grad </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cr)</span></span></code></pre></div><p><img src="'+E+`" alt=""></p><h2 id="Plot-the-effect-of-the-new-liquid-kr-exponent-on-the-gas-production" tabindex="-1">Plot the effect of the new liquid kr exponent on the gas production <a class="header-anchor" href="#Plot-the-effect-of-the-new-liquid-kr-exponent-on-the-gas-production" aria-label="Permalink to &quot;Plot the effect of the new liquid kr exponent on the gas production {#Plot-the-effect-of-the-new-liquid-kr-exponent-on-the-gas-production}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> big</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-7</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-7</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">else</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">8e-6</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">8e-6</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">kre </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data_domain_with_gradients[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:KrExponents</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">exp_l </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> kre[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;exp_liquid&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, exp_l, is_grad </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cr)</span></span></code></pre></div><p><img src="`+c+`" alt=""></p><h2 id="Plot-the-effect-of-the-new-vapor-kr-exponent-on-the-gas-production" tabindex="-1">Plot the effect of the new vapor kr exponent on the gas production <a class="header-anchor" href="#Plot-the-effect-of-the-new-vapor-kr-exponent-on-the-gas-production" aria-label="Permalink to &quot;Plot the effect of the new vapor kr exponent on the gas production {#Plot-the-effect-of-the-new-vapor-kr-exponent-on-the-gas-production}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">exp_v </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> kre[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;exp_vapor&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, exp_v, is_grad </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cr)</span></span></code></pre></div><p><img src="`+d+`" alt=""></p><h2 id="Plot-the-effect-of-the-liquid-phase-viscosity" tabindex="-1">Plot the effect of the liquid phase viscosity <a class="header-anchor" href="#Plot-the-effect-of-the-liquid-phase-viscosity" aria-label="Permalink to &quot;Plot the effect of the liquid phase viscosity {#Plot-the-effect-of-the-liquid-phase-viscosity}&quot;">​</a></h2><p>Note: The viscosity can in many models be a variable and not a parameter. For this simple model, however, it is treated as a parameter and we obtain sensitivities.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">mu </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data_domain_with_gradients[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:PhaseViscosities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> big</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.001</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">else</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.01</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.01</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">mu_l </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> mu[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;mu_liquid&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, mu_l, is_grad </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cr)</span></span></code></pre></div><p><img src="`+g+`" alt=""></p><h2 id="Plot-the-effect-of-the-liquid-phase-viscosity-2" tabindex="-1">Plot the effect of the liquid phase viscosity <a class="header-anchor" href="#Plot-the-effect-of-the-liquid-phase-viscosity-2" aria-label="Permalink to &quot;Plot the effect of the liquid phase viscosity {#Plot-the-effect-of-the-liquid-phase-viscosity-2}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">mu_v </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> mu[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">myplot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;mu_vapor&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, mu_v, is_grad </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cr)</span></span></code></pre></div><p><img src="`+y+'" alt=""></p><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/intro_sensitivities.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/intro_sensitivities.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',63)]))}const f=a(o,[["render",u]]);export{D as __pageData,f as default};
