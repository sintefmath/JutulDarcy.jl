import{_ as a,c as i,a5 as n,o as l}from"./chunks/framework.BnJvuozR.js";const p="/JutulDarcy.jl/previews/PR65/assets/pjouopj.eZtJR1xw.jpeg",e="/JutulDarcy.jl/previews/PR65/assets/cwuughl.B7Is6nt6.jpeg",t="/JutulDarcy.jl/previews/PR65/assets/elzqayp.D49xvrLY.jpeg",h="/JutulDarcy.jl/previews/PR65/assets/kehdvhl.0FcmmgE8.jpeg",k="/JutulDarcy.jl/previews/PR65/assets/llazgrc.DMhZgfRz.jpeg",r="/JutulDarcy.jl/previews/PR65/assets/vkjkvmn.DxEcouAp.jpeg",o="/JutulDarcy.jl/previews/PR65/assets/zagpzpw.B3zQ-Tkt.jpeg",E="/JutulDarcy.jl/previews/PR65/assets/ntquwwy.CpaXmJvq.jpeg",c="/JutulDarcy.jl/previews/PR65/assets/gigpuqj.Bjq_6Nj2.jpeg",D=JSON.parse('{"title":"Norne: Real field black-oil model","description":"","frontmatter":{},"headers":[],"relativePath":"examples/validation_norne_nohyst.md","filePath":"examples/validation_norne_nohyst.md","lastUpdated":null}'),d={name:"examples/validation_norne_nohyst.md"};function y(g,s,F,u,m,C){return l(),i("div",null,s[0]||(s[0]=[n(`<h1 id="Norne:-Real-field-black-oil-model" tabindex="-1">Norne: Real field black-oil model <a class="header-anchor" href="#Norne:-Real-field-black-oil-model" aria-label="Permalink to &quot;Norne: Real field black-oil model {#Norne:-Real-field-black-oil-model}&quot;">​</a></h1><p>The Norne model is a real field model. The model has been adapted so that the input file only contains features present in JutulDarcy, with the most notable omissions being removal of hysteresis and threshold pressures between equilibriation reqgions. For more details, see the <a href="https://opm-project.org/?page_id=559" target="_blank" rel="noreferrer">OPM data webpage</a></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Jutul, JutulDarcy, GLMakie, DelimitedFiles, HYPRE, GeoEnergyIO</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">norne_dir </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeoEnergyIO</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">test_input_file_path</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;NORNE_NOHYST&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">data_pth </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> joinpath</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(norne_dir, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;NORNE_NOHYST.DATA&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">data </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> parse_data_file</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(data_pth)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">case </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_case_from_data_file</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(data)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Jutul case with 247 time-steps (9 years, 3 weeks, 3.817 days) and forces for each step.</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Model:</span></span>
<span class="line"><span></span></span>
<span class="line"><span>MultiModel with 38 models and 108 cross-terms. 133395 equations, 133395 degrees of freedom and 1068242 parameters.</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  models:</span></span>
<span class="line"><span>    1) Reservoir (133251x133251)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ MinimalTPFATopology (44417 cells, 133847 faces)</span></span>
<span class="line"><span>    2) B-4AH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [B-4AH] (1 nodes, 0 segments, 29 perforations)</span></span>
<span class="line"><span>    3) B-2H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [B-2H] (1 nodes, 0 segments, 8 perforations)</span></span>
<span class="line"><span>    4) E-3H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [E-3H] (1 nodes, 0 segments, 21 perforations)</span></span>
<span class="line"><span>    5) D-3BH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [D-3BH] (1 nodes, 0 segments, 13 perforations)</span></span>
<span class="line"><span>    6) E-3BH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [E-3BH] (1 nodes, 0 segments, 26 perforations)</span></span>
<span class="line"><span>    7) B-1BH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [B-1BH] (1 nodes, 0 segments, 6 perforations)</span></span>
<span class="line"><span>    8) D-4AH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [D-4AH] (1 nodes, 0 segments, 6 perforations)</span></span>
<span class="line"><span>    9) B-4BH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [B-4BH] (1 nodes, 0 segments, 10 perforations)</span></span>
<span class="line"><span>    10) C-4AH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [C-4AH] (1 nodes, 0 segments, 15 perforations)</span></span>
<span class="line"><span>    11) E-3CH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [E-3CH] (1 nodes, 0 segments, 8 perforations)</span></span>
<span class="line"><span>    12) B-1AH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [B-1AH] (1 nodes, 0 segments, 5 perforations)</span></span>
<span class="line"><span>    13) D-4H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [D-4H] (1 nodes, 0 segments, 26 perforations)</span></span>
<span class="line"><span>    14) E-1H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [E-1H] (1 nodes, 0 segments, 10 perforations)</span></span>
<span class="line"><span>    15) B-4H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [B-4H] (1 nodes, 0 segments, 20 perforations)</span></span>
<span class="line"><span>    16) E-4H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [E-4H] (1 nodes, 0 segments, 10 perforations)</span></span>
<span class="line"><span>    17) K-3H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [K-3H] (1 nodes, 0 segments, 11 perforations)</span></span>
<span class="line"><span>    18) C-2H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [C-2H] (1 nodes, 0 segments, 8 perforations)</span></span>
<span class="line"><span>    19) F-1H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [F-1H] (1 nodes, 0 segments, 18 perforations)</span></span>
<span class="line"><span>    20) E-2H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [E-2H] (1 nodes, 0 segments, 9 perforations)</span></span>
<span class="line"><span>    21) B-3H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [B-3H] (1 nodes, 0 segments, 14 perforations)</span></span>
<span class="line"><span>    22) F-4H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [F-4H] (1 nodes, 0 segments, 17 perforations)</span></span>
<span class="line"><span>    23) D-1H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [D-1H] (1 nodes, 0 segments, 8 perforations)</span></span>
<span class="line"><span>    24) F-3H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [F-3H] (1 nodes, 0 segments, 23 perforations)</span></span>
<span class="line"><span>    25) F-2H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [F-2H] (1 nodes, 0 segments, 20 perforations)</span></span>
<span class="line"><span>    26) C-3H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [C-3H] (1 nodes, 0 segments, 21 perforations)</span></span>
<span class="line"><span>    27) D-2H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [D-2H] (1 nodes, 0 segments, 9 perforations)</span></span>
<span class="line"><span>    28) B-1H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [B-1H] (1 nodes, 0 segments, 16 perforations)</span></span>
<span class="line"><span>    29) D-3H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [D-3H] (1 nodes, 0 segments, 24 perforations)</span></span>
<span class="line"><span>    30) D-3AH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [D-3AH] (1 nodes, 0 segments, 12 perforations)</span></span>
<span class="line"><span>    31) D-1CH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [D-1CH] (1 nodes, 0 segments, 15 perforations)</span></span>
<span class="line"><span>    32) B-4DH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [B-4DH] (1 nodes, 0 segments, 11 perforations)</span></span>
<span class="line"><span>    33) E-4AH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [E-4AH] (1 nodes, 0 segments, 12 perforations)</span></span>
<span class="line"><span>    34) C-4H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [C-4H] (1 nodes, 0 segments, 7 perforations)</span></span>
<span class="line"><span>    35) E-3AH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [E-3AH] (1 nodes, 0 segments, 7 perforations)</span></span>
<span class="line"><span>    36) C-1H (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [C-1H] (1 nodes, 0 segments, 21 perforations)</span></span>
<span class="line"><span>    37) E-2AH (3x3)</span></span>
<span class="line"><span>       StandardBlackOilSystem with (AqueousPhase(), LiquidPhase(), VaporPhase())</span></span>
<span class="line"><span>       ∈ SimpleWell [E-2AH] (1 nodes, 0 segments, 7 perforations)</span></span>
<span class="line"><span>    38) Facility (36x36)</span></span>
<span class="line"><span>       JutulDarcy.PredictionMode()</span></span>
<span class="line"><span>       ∈ WellGroup([Symbol(&quot;B-4AH&quot;), Symbol(&quot;B-2H&quot;), Symbol(&quot;E-3H&quot;), Symbol(&quot;D-3BH&quot;), Symbol(&quot;E-3BH&quot;), Symbol(&quot;B-1BH&quot;), Symbol(&quot;D-4AH&quot;), Symbol(&quot;B-4BH&quot;), Symbol(&quot;C-4AH&quot;), Symbol(&quot;E-3CH&quot;), Symbol(&quot;B-1AH&quot;), Symbol(&quot;D-4H&quot;), Symbol(&quot;E-1H&quot;), Symbol(&quot;B-4H&quot;), Symbol(&quot;E-4H&quot;), Symbol(&quot;K-3H&quot;), Symbol(&quot;C-2H&quot;), Symbol(&quot;F-1H&quot;), Symbol(&quot;E-2H&quot;), Symbol(&quot;B-3H&quot;), Symbol(&quot;F-4H&quot;), Symbol(&quot;D-1H&quot;), Symbol(&quot;F-3H&quot;), Symbol(&quot;F-2H&quot;), Symbol(&quot;C-3H&quot;), Symbol(&quot;D-2H&quot;), Symbol(&quot;B-1H&quot;), Symbol(&quot;D-3H&quot;), Symbol(&quot;D-3AH&quot;), Symbol(&quot;D-1CH&quot;), Symbol(&quot;B-4DH&quot;), Symbol(&quot;E-4AH&quot;), Symbol(&quot;C-4H&quot;), Symbol(&quot;E-3AH&quot;), Symbol(&quot;C-1H&quot;), Symbol(&quot;E-2AH&quot;)], true, true)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  cross_terms:</span></span>
<span class="line"><span>    1) B-4AH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    2) B-2H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    3) E-3H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    4) D-3BH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    5) E-3BH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    6) B-1BH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    7) D-4AH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    8) B-4BH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    9) C-4AH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    10) E-3CH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    11) B-1AH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    12) D-4H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    13) E-1H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    14) B-4H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    15) E-4H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    16) K-3H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    17) C-2H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    18) F-1H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    19) E-2H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    20) B-3H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    21) F-4H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    22) D-1H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    23) F-3H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    24) F-2H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    25) C-3H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    26) D-2H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    27) B-1H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    28) D-3H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    29) D-3AH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    30) D-1CH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    31) B-4DH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    32) E-4AH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    33) C-4H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    34) E-3AH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    35) C-1H &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    36) E-2AH &lt;-&gt; Reservoir (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.ReservoirFromWellFlowCT</span></span>
<span class="line"><span>    37) B-4AH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    38) Facility  -&gt; B-4AH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    39) B-2H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    40) Facility  -&gt; B-2H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    41) E-3H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    42) Facility  -&gt; E-3H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    43) D-3BH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    44) Facility  -&gt; D-3BH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    45) E-3BH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    46) Facility  -&gt; E-3BH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    47) B-1BH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    48) Facility  -&gt; B-1BH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    49) D-4AH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    50) Facility  -&gt; D-4AH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    51) B-4BH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    52) Facility  -&gt; B-4BH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    53) C-4AH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    54) Facility  -&gt; C-4AH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    55) E-3CH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    56) Facility  -&gt; E-3CH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    57) B-1AH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    58) Facility  -&gt; B-1AH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    59) D-4H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    60) Facility  -&gt; D-4H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    61) E-1H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    62) Facility  -&gt; E-1H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    63) B-4H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    64) Facility  -&gt; B-4H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    65) E-4H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    66) Facility  -&gt; E-4H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    67) K-3H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    68) Facility  -&gt; K-3H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    69) C-2H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    70) Facility  -&gt; C-2H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    71) F-1H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    72) Facility  -&gt; F-1H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    73) E-2H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    74) Facility  -&gt; E-2H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    75) B-3H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    76) Facility  -&gt; B-3H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    77) F-4H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    78) Facility  -&gt; F-4H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    79) D-1H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    80) Facility  -&gt; D-1H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    81) F-3H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    82) Facility  -&gt; F-3H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    83) F-2H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    84) Facility  -&gt; F-2H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    85) C-3H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    86) Facility  -&gt; C-3H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    87) D-2H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    88) Facility  -&gt; D-2H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    89) B-1H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    90) Facility  -&gt; B-1H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    91) D-3H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    92) Facility  -&gt; D-3H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    93) D-3AH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    94) Facility  -&gt; D-3AH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    95) D-1CH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    96) Facility  -&gt; D-1CH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    97) B-4DH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    98) Facility  -&gt; B-4DH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    99) E-4AH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    100) Facility  -&gt; E-4AH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    101) C-4H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    102) Facility  -&gt; C-4H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    103) E-3AH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    104) Facility  -&gt; E-3AH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    105) C-1H  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    106) Facility  -&gt; C-1H (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span>    107) E-2AH  -&gt; Facility (Eq: control_equation)</span></span>
<span class="line"><span>       JutulDarcy.FacilityFromWellFlowCT</span></span>
<span class="line"><span>    108) Facility  -&gt; E-2AH (Eq: mass_conservation)</span></span>
<span class="line"><span>       JutulDarcy.WellFromFacilityFlowCT</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Model storage will be optimized for runtime performance.</span></span></code></pre></div><h2 id="Unpack-the-case-to-see-basic-data-structures" tabindex="-1">Unpack the case to see basic data structures <a class="header-anchor" href="#Unpack-the-case-to-see-basic-data-structures" aria-label="Permalink to &quot;Unpack the case to see basic data structures {#Unpack-the-case-to-see-basic-data-structures}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">model </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> case</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">model</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">parameters </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> case</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">parameters</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">forces </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> case</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">forces</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> case</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dt;</span></span></code></pre></div><h2 id="Plot-the-reservoir-mesh,-wells-and-faults" tabindex="-1">Plot the reservoir mesh, wells and faults <a class="header-anchor" href="#Plot-the-reservoir-mesh,-wells-and-faults" aria-label="Permalink to &quot;Plot the reservoir mesh, wells and faults {#Plot-the-reservoir-mesh,-wells-and-faults}&quot;">​</a></h2><p>We compose a few different plotting calls together to make a plot that shows the outline of the mesh, the fault structures and the well trajectories.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">import</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Jutul</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> plot_mesh_edges!</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">reservoir </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reservoir_domain</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">mesh </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> physical_representation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(reservoir)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">wells </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> get_model_wells</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1200</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">800</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], zreversed </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_mesh_edges!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, mesh, alpha </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (k, w) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> wells</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    plot_well!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, mesh, w)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_faults!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, mesh, alpha </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">azimuth[] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3.0</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">elevation[] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.5</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+p+`" alt=""></p><h2 id="Plot-the-reservoir-static-properties-in-interactive-viewer" tabindex="-1">Plot the reservoir static properties in interactive viewer <a class="header-anchor" href="#Plot-the-reservoir-static-properties-in-interactive-viewer" aria-label="Permalink to &quot;Plot the reservoir static properties in interactive viewer {#Plot-the-reservoir-static-properties-in-interactive-viewer}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> plot_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, key </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :porosity</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> fig</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">current_axis[]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_faults!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, mesh, alpha </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">azimuth[] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3.0</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">elevation[] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.5</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+e+`" alt=""></p><h2 id="Simulate-the-model" tabindex="-1">Simulate the model <a class="header-anchor" href="#Simulate-the-model" aria-label="Permalink to &quot;Simulate the model {#Simulate-the-model}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ws, states </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case, output_substates </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>ReservoirSimResult with 791 entries:</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  wells (36 present):</span></span>
<span class="line"><span>    :B-4H</span></span>
<span class="line"><span>    :F-4H</span></span>
<span class="line"><span>    :D-4H</span></span>
<span class="line"><span>    :F-2H</span></span>
<span class="line"><span>    :D-2H</span></span>
<span class="line"><span>    :B-1H</span></span>
<span class="line"><span>    :C-4H</span></span>
<span class="line"><span>    :C-1H</span></span>
<span class="line"><span>    :B-2H</span></span>
<span class="line"><span>    :E-1H</span></span>
<span class="line"><span>    :B-4BH</span></span>
<span class="line"><span>    :D-3AH</span></span>
<span class="line"><span>    :D-3BH</span></span>
<span class="line"><span>    :C-3H</span></span>
<span class="line"><span>    :E-3H</span></span>
<span class="line"><span>    :K-3H</span></span>
<span class="line"><span>    :E-4AH</span></span>
<span class="line"><span>    :D-1H</span></span>
<span class="line"><span>    :B-1AH</span></span>
<span class="line"><span>    :E-3CH</span></span>
<span class="line"><span>    :E-4H</span></span>
<span class="line"><span>    :E-2H</span></span>
<span class="line"><span>    :E-2AH</span></span>
<span class="line"><span>    :C-4AH</span></span>
<span class="line"><span>    :B-1BH</span></span>
<span class="line"><span>    :C-2H</span></span>
<span class="line"><span>    :B-4DH</span></span>
<span class="line"><span>    :D-3H</span></span>
<span class="line"><span>    :E-3AH</span></span>
<span class="line"><span>    :D-4AH</span></span>
<span class="line"><span>    :B-3H</span></span>
<span class="line"><span>    :F-3H</span></span>
<span class="line"><span>    :E-3BH</span></span>
<span class="line"><span>    :D-1CH</span></span>
<span class="line"><span>    :F-1H</span></span>
<span class="line"><span>    :B-4AH</span></span>
<span class="line"><span>    Results per well:</span></span>
<span class="line"><span>       :wrat =&gt; Vector{Float64} of size (791,)</span></span>
<span class="line"><span>       :Aqueous_mass_rate =&gt; Vector{Float64} of size (791,)</span></span>
<span class="line"><span>       :orat =&gt; Vector{Float64} of size (791,)</span></span>
<span class="line"><span>       :bhp =&gt; Vector{Float64} of size (791,)</span></span>
<span class="line"><span>       :lrat =&gt; Vector{Float64} of size (791,)</span></span>
<span class="line"><span>       :mass_rate =&gt; Vector{Float64} of size (791,)</span></span>
<span class="line"><span>       :rate =&gt; Vector{Float64} of size (791,)</span></span>
<span class="line"><span>       :Vapor_mass_rate =&gt; Vector{Float64} of size (791,)</span></span>
<span class="line"><span>       :control =&gt; Vector{Symbol} of size (791,)</span></span>
<span class="line"><span>       :Liquid_mass_rate =&gt; Vector{Float64} of size (791,)</span></span>
<span class="line"><span>       :grat =&gt; Vector{Float64} of size (791,)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  states (Vector with 791 entries, reservoir variables for each state)</span></span>
<span class="line"><span>    :Rv =&gt; Vector{Float64} of size (44417,)</span></span>
<span class="line"><span>    :BlackOilUnknown =&gt; Vector{BlackOilX{Float64}} of size (44417,)</span></span>
<span class="line"><span>    :Saturations =&gt; Matrix{Float64} of size (3, 44417)</span></span>
<span class="line"><span>    :Pressure =&gt; Vector{Float64} of size (44417,)</span></span>
<span class="line"><span>    :Rs =&gt; Vector{Float64} of size (44417,)</span></span>
<span class="line"><span>    :ImmiscibleSaturation =&gt; Vector{Float64} of size (44417,)</span></span>
<span class="line"><span>    :TotalMasses =&gt; Matrix{Float64} of size (3, 44417)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  time (report time for each state)</span></span>
<span class="line"><span>     Vector{Float64} of length 791</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  result (extended states, reports)</span></span>
<span class="line"><span>     SimResult with 247 entries</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  extra</span></span>
<span class="line"><span>     Dict{Any, Any} with keys :simulator, :config</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Completed at Sep. 27 2024 20:09 after 22 minutes, 59 seconds, 584.2 milliseconds.</span></span></code></pre></div><h2 id="Plot-the-reservoir-solution" tabindex="-1">Plot the reservoir solution <a class="header-anchor" href="#Plot-the-reservoir-solution" aria-label="Permalink to &quot;Plot the reservoir solution {#Plot-the-reservoir-solution}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> plot_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, states, step </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 247</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, key </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Saturations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> fig</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">current_axis[]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">azimuth[] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3.0</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">elevation[] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.5</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+t+`" alt=""></p><h2 id="Load-reference-and-set-up-plotting" tabindex="-1">Load reference and set up plotting <a class="header-anchor" href="#Load-reference-and-set-up-plotting" aria-label="Permalink to &quot;Load reference and set up plotting {#Load-reference-and-set-up-plotting}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">csv_path </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> joinpath</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(norne_dir, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;REFERENCE.CSV&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">data, header </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> readdlm</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(csv_path, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&#39;,&#39;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, header </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">time_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data[:, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">time_jutul </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> deepcopy</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ws</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">time)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">wells </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> deepcopy</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ws</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">wells)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">wnames </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> collect</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">keys</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(wells))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">nw </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> length</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(wnames)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">day </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> si_unit</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:day</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">cmap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :tableau_hue_circle</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">inj </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Symbol[]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">prod </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Symbol[]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (wellname, well) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> pairs</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(wells)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    qts </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> well[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:wrat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> well[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:orat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> well[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:grat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sum</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(qts) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        push!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(inj, wellname)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    else</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        push!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prod, wellname)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> plot_well_comparison</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(response, well_names, reponse_name </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$response</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1000</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">400</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> response </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :bhp</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        ys </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">si_unit</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:bar</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        yl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Bottom hole pressure / Bar&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    elseif</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> response </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :wrat</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        ys </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> si_unit</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:day</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        yl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Surface water rate / m³/day&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    elseif</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> response </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :grat</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        ys </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> si_unit</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:day</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e6</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        yl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Surface gas rate / 10⁶ m³/day&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    elseif</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> response </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :orat</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        ys </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> si_unit</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:day</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">/</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1000</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">si_unit</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:stb</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        yl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Surface oil rate / 10³ stb/day&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    else</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        error</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$response</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> not ready.&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    welltypes </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> []</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], xlabel </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Time / days&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, ylabel </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> yl)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    linehandles </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> []</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    linelabels </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> []</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> well_name </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> well_names</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        well </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> wells[well_name]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        label_in_csv </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$well_name</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$response</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        ref_pos </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> findfirst</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> label_in_csv, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">vec</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(header))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        qoi </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> copy</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(well[response])</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ys</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        qoi_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data[:, ref_pos]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ys</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        tot_rate </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> copy</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(well[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:rate</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        @.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> qoi[tot_rate </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> NaN</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        grat_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data[:, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">findfirst</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$well_name</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">:grat&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">vec</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(header))]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        orat_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data[:, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">findfirst</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$well_name</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">:orat&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">vec</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(header))]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        wrat_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data[:, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">findfirst</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$well_name</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">:wrat&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">vec</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(header))]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        tot_rate_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> grat_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> orat_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> wrat_ref</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        @.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> qoi_ref[tot_rate_ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> NaN</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        crange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">max</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">length</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(well_names), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        lh </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, time_jutul</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">./</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">day, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">abs</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(qoi),</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> crange,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            label </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$well_name</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, colormap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cmap</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        )</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        push!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(linehandles, lh)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        push!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(linelabels, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$well_name</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, time_ref</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">./</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">day, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">abs</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(qoi_ref),</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> crange,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            linestyle </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :dash</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            colormap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cmap</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        )</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    l1 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> LineElement</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :black</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, linestyle </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> nothing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    l2 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> LineElement</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :black</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, linestyle </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :dash</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    Legend</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], linehandles, linelabels, nbanks </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    Legend</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], [l1, l2], [</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;JutulDarcy.jl&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;OPM Flow&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    fig</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot_well_comparison (generic function with 2 methods)</span></span></code></pre></div><h2 id="Injector-bhp" tabindex="-1">Injector bhp <a class="header-anchor" href="#Injector-bhp" aria-label="Permalink to &quot;Injector bhp {#Injector-bhp}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_well_comparison</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:bhp</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, inj, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Bottom hole pressure&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="`+h+'" alt=""></p><h2 id="Gas-injection-rates" tabindex="-1">Gas injection rates <a class="header-anchor" href="#Gas-injection-rates" aria-label="Permalink to &quot;Gas injection rates {#Gas-injection-rates}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_well_comparison</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:grat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, inj, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Gas surface injection rate&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="'+k+'" alt=""></p><h2 id="Water-injection-rates" tabindex="-1">Water injection rates <a class="header-anchor" href="#Water-injection-rates" aria-label="Permalink to &quot;Water injection rates {#Water-injection-rates}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_well_comparison</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:wrat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, inj, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Water surface injection rate&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="'+r+'" alt=""></p><h2 id="Producer-bhp" tabindex="-1">Producer bhp <a class="header-anchor" href="#Producer-bhp" aria-label="Permalink to &quot;Producer bhp {#Producer-bhp}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_well_comparison</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:bhp</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, prod, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Bottom hole pressure&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="'+o+'" alt=""></p><h2 id="Oil-production-rates" tabindex="-1">Oil production rates <a class="header-anchor" href="#Oil-production-rates" aria-label="Permalink to &quot;Oil production rates {#Oil-production-rates}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_well_comparison</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:orat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, prod, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Oil surface production rate&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="'+E+'" alt=""></p><h2 id="Gas-production-rates" tabindex="-1">Gas production rates <a class="header-anchor" href="#Gas-production-rates" aria-label="Permalink to &quot;Gas production rates {#Gas-production-rates}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_well_comparison</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:grat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, prod, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Gas surface production rate&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="'+c+'" alt=""></p><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/validation_norne_nohyst.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/validation_norne_nohyst.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',44)]))}const B=a(d,[["render",y]]);export{D as __pageData,B as default};
