import{_ as n,c as r,a5 as a,j as s,a as t,G as l,B as o,o as p}from"./chunks/framework.CJGehonb.js";const B=JSON.parse('{"title":"Utilities","description":"","frontmatter":{},"headers":[],"relativePath":"man/basics/utilities.md","filePath":"man/basics/utilities.md","lastUpdated":null}'),h={name:"man/basics/utilities.md"},d={class:"jldocstring custom-block",open:""},u={class:"jldocstring custom-block",open:""},c={class:"jldocstring custom-block",open:""},k={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""},y={class:"jldocstring custom-block",open:""},b={class:"jldocstring custom-block",open:""},m={class:"jldocstring custom-block",open:""},f={class:"jldocstring custom-block",open:""},E={class:"jldocstring custom-block",open:""},v={class:"jldocstring custom-block",open:""},C={class:"jldocstring custom-block",open:""},F={class:"jldocstring custom-block",open:""},j={class:"jldocstring custom-block",open:""},_={class:"jldocstring custom-block",open:""};function D(O,i,w,x,P,J){const e=o("Badge");return p(),r("div",null,[i[45]||(i[45]=a('<h1 id="utilities" tabindex="-1">Utilities <a class="header-anchor" href="#utilities" aria-label="Permalink to &quot;Utilities&quot;">​</a></h1><p>This section describes various utilities that do not fit in other sections.</p><h2 id="CO2-and-brine-correlations" tabindex="-1">CO2 and brine correlations <a class="header-anchor" href="#CO2-and-brine-correlations" aria-label="Permalink to &quot;CO2 and brine correlations {#CO2-and-brine-correlations}&quot;">​</a></h2><p>These functions are not exported, but can be found inside the <code>CO2Properties</code> submodule. The functions described here are a Julia port of the MRST module described in [<a href="/JutulDarcy.jl/v0.2.35/extras/refs#salo_co2">6</a>] They can be accessed by explicit import:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">import</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CO2Properties</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> name_of_function</span></span></code></pre></div>',5)),s("details",d,[s("summary",null,[i[0]||(i[0]=s("a",{id:"JutulDarcy.CO2Properties.compute_co2_brine_props",href:"#JutulDarcy.CO2Properties.compute_co2_brine_props"},[s("span",{class:"jlbinding"},"JutulDarcy.CO2Properties.compute_co2_brine_props")],-1)),i[1]||(i[1]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[2]||(i[2]=a(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">props </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> compute_co2_brine_props</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(p_pascal, T_K, salt_mole_fractions </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Float64[], salt_names </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> String[];</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    check</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    iterate</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    maxits</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">15</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ionized</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">false</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">props </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> compute_co2_brine_props</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">200e5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">273.15</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> +</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 30.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Get pure phase properties (density, viscosity) and equilibrium constants for brine and water. This functions and the functions used in this function are heavily based on comparable code in the the co2lab-mit MRST module developed by Lluis Salo.</p><p>The salt mole fractions are of the total brine (i.e. including H2O in the calculation) and can any subset of the following names provided:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>(&quot;NaCl&quot;, &quot;KCl&quot;, &quot;CaSO4&quot;, &quot;CaCl2&quot;, &quot;MgSO4&quot;, &quot;MgCl2&quot;)</span></span></code></pre></div><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/CO2Properties/generation.jl#L528-L547" target="_blank" rel="noreferrer">source</a></p>`,5))]),s("details",u,[s("summary",null,[i[3]||(i[3]=s("a",{id:"JutulDarcy.CO2Properties.pvt_brine_RoweChou1970",href:"#JutulDarcy.CO2Properties.pvt_brine_RoweChou1970"},[s("span",{class:"jlbinding"},"JutulDarcy.CO2Properties.pvt_brine_RoweChou1970")],-1)),i[4]||(i[4]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[5]||(i[5]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[rho_brine, c_brine] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> pvt_brine_RoweChou1970</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(T, P, S)</span></span></code></pre></div><p>Calculate brine density and/or compressibility using Rowe and Chou&quot;s (1970) correlation.</p><p><strong>Parameter range for correlation:</strong></p><p>P &lt;= 35 MPa 293.15 &lt;= T &lt;= 423.15 K</p><p><strong>Arguments:</strong></p><ul><li><p>T: Scalar with temperature value in Kelvin</p></li><li><p>P: Scalar with pressure value in bar</p></li><li><p>S: Salt mass fraction</p></li></ul><p><strong>Outputs</strong></p><ul><li><p>rho_brine: Scalar with density value in kg/m3</p></li><li><p>c_brine: Scalar with compressibility value in 1/kPa</p></li></ul><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/CO2Properties/generation.jl#L1-L23" target="_blank" rel="noreferrer">source</a></p>',9))]),s("details",c,[s("summary",null,[i[6]||(i[6]=s("a",{id:"JutulDarcy.CO2Properties.activity_co2_DS2003",href:"#JutulDarcy.CO2Properties.activity_co2_DS2003"},[s("span",{class:"jlbinding"},"JutulDarcy.CO2Properties.activity_co2_DS2003")],-1)),i[7]||(i[7]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[8]||(i[8]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">gamma_co2 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> activity_co2_DS2003</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(T, P, m_io)</span></span></code></pre></div><p>Calculate a CO2 pseudo activity coefficient based on a virial expansion of excess Gibbs energy.</p><p><strong>Arguments</strong></p><ul><li><p>T: Scalar with temperature value in Kelvin</p></li><li><p>P: Scalar with pressure value in bar</p></li><li><p>m_io: Vector where each entry corresponds to the molality of a particular ion in the initial brine solution. The order is as follows: [ Na(+), K(+), Ca(2+), Mg(2+), Cl(-), SO4(2-)]</p></li></ul><p><strong>Outputs</strong></p><ul><li><p>V_m: Scalar with molar volume in [cm^3/mol]</p></li><li><p>rhox: Scalar with density in [mol/m^3]</p></li><li><p>rho: Scalar with density in [kg/m^3]</p></li></ul><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/CO2Properties/generation.jl#L263-L281" target="_blank" rel="noreferrer">source</a></p>',7))]),s("details",k,[s("summary",null,[i[9]||(i[9]=s("a",{id:"JutulDarcy.CO2Properties.viscosity_brine_co2_mixture_IC2012",href:"#JutulDarcy.CO2Properties.viscosity_brine_co2_mixture_IC2012"},[s("span",{class:"jlbinding"},"JutulDarcy.CO2Properties.viscosity_brine_co2_mixture_IC2012")],-1)),i[10]||(i[10]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[11]||(i[11]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">mu_b_co2 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> viscosity_brine_co2_mixture_IC2012</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(T, P, m_nacl, w_co2)</span></span></code></pre></div><p>Calculate the dynamic viscosity of a solution of H2O + NaCl (brine) with dissolved CO2.</p><p><strong>Parameter range for correlation:</strong></p><p>For pure water + CO2, the model is based on experimental data by Kumagai et al. (1998), which is for p up to 400 bar and T up to 50 C, and Bando et al. (2004), which is valid in 30-60C and 10-20MPa. The model of Mao &amp; Duan (2009) for brine viscosity reaches 623K, 1000 bar and high ionic strength. However, the model used to determine the viscosity when co2 dissolves in the brine (Islam &amp; Carlson, 2012) is based on experimental data by Bando et al. (2004) and Fleury and Deschamps (2008), who provided experimental data up to P = 200 bar, T extrapolated to 100 C, and maximum salinity of 2.7M.</p><p><strong>Arguments</strong></p><ul><li><p>T: Scalar with temperature value in Kelvin</p></li><li><p>P: Scalar with pressure value in bar</p></li><li><p>m_nacl: Salt molality (NaCl) in mol/kg solvent</p></li><li><p>w_co2: Mass fraction of CO2 in the aqueous solution (i.e. brine)</p></li></ul><p><strong>Outputs</strong></p><ul><li>mu_b_co2: Scalar with dynamic viscosity in Pa*s</li></ul><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/CO2Properties/generation.jl#L365-L394" target="_blank" rel="noreferrer">source</a></p>',9))]),s("details",g,[s("summary",null,[i[12]||(i[12]=s("a",{id:"JutulDarcy.CO2Properties.pvt_brine_BatzleWang1992",href:"#JutulDarcy.CO2Properties.pvt_brine_BatzleWang1992"},[s("span",{class:"jlbinding"},"JutulDarcy.CO2Properties.pvt_brine_BatzleWang1992")],-1)),i[13]||(i[13]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[14]||(i[14]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rho_b </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> pvt_brine_BatzleWang1992</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(T, P, w_nacl)</span></span></code></pre></div><p>Calculate the brine (H2O + NaCl) density based on Batzle &amp; Wang (1992). These authors used the data of Rowe and Chou (1970), Zarembo &amp; Fedorov (1975) and Potter &amp; Brown (1977) to expand the P, T validity range.</p><p><strong>Parameter range for correlation:</strong></p><p>P valid from 5 to 100 MPa, T from 20 to 350 C (Adams &amp; Bachu, 2002)</p><p><strong>Arguments</strong></p><ul><li><p>T: Temperature value in Kelvin</p></li><li><p>P: Pressure value in bar</p></li><li><p>w_nacl: Salt (NaCl) mass fraction</p></li></ul><p><strong>Outputs</strong></p><ul><li>rho_b: Scalar with brine density in kg/m3</li></ul><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/CO2Properties/generation.jl#L110-L129" target="_blank" rel="noreferrer">source</a></p>',9))]),s("details",y,[s("summary",null,[i[15]||(i[15]=s("a",{id:"JutulDarcy.CO2Properties.viscosity_co2_Fenghour1998",href:"#JutulDarcy.CO2Properties.viscosity_co2_Fenghour1998"},[s("span",{class:"jlbinding"},"JutulDarcy.CO2Properties.viscosity_co2_Fenghour1998")],-1)),i[16]||(i[16]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[17]||(i[17]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">mu </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> viscosity_co2_Fenghour1998</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(T, rho)</span></span></code></pre></div><p>Calculate CO2 viscosity from Vesovic et al., J Phys Chem Ref Data (1990) and Fenghour et al., J Phys Chem Ref Data (1998), as described in Hassanzadeh et al., IJGGC (2008).</p><p><strong>Arguments:</strong></p><ul><li><p>T: Scalar of temperature value in Kelvin</p></li><li><p>rho: Scalar of density value in kg/m^3</p></li></ul><p><strong>Outputs</strong></p><ul><li>mu: Dynamic viscosity in Pa*s</li></ul><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/CO2Properties/generation.jl#L60-L76" target="_blank" rel="noreferrer">source</a></p>',7))]),s("details",b,[s("summary",null,[i[18]||(i[18]=s("a",{id:"JutulDarcy.CO2Properties.pvt_co2_RedlichKwong1949",href:"#JutulDarcy.CO2Properties.pvt_co2_RedlichKwong1949"},[s("span",{class:"jlbinding"},"JutulDarcy.CO2Properties.pvt_co2_RedlichKwong1949")],-1)),i[19]||(i[19]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[20]||(i[20]=a(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[V_m, rhox, rho] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> pvt_co2_RedlichKwong1949</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(T, P)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[V_m, rhox, rho] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> pvt_co2_RedlichKwong1949</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(T, P, a_m, b_m)</span></span></code></pre></div><p>Calculate CO2 molar volume and density using Redlich and Kwong (1949) EoS (= RK EoS).</p><p><strong>Parameter range for correlation:</strong></p><p>Tested by Spycher et al. (2003) with constant intermolecular attraction parameters to yield accurate results in the T range ~10 to ~100C and P range up to 600 bar, for (1) CO2 compressibility factor, (2) CO2 fugacity coefficient and (3) mutual solubilities of H2O and CO2 in the gas and aqueous phase (respectively).</p><p><strong>Arguments</strong></p><ul><li><p>T: Scalar with temperature value in Kelvin</p></li><li><p>P: Scalar with pressure value in bar</p></li></ul><p><strong>Optional arguments:</strong></p><ul><li><p>a_m: Intermolecular attraction constant (of the mixture) in bar_cm^6_K^0.5/mol^2</p></li><li><p>b_m: Intermolecular repulsion constant (of the mixture) in cm^3/mol</p></li></ul><p><strong>Outputs</strong></p><ul><li><p>V_m: Scalar with molar volume in [cm^3/mol]</p></li><li><p>rhox: Scalar with density in [mol/m^3]</p></li><li><p>rho: Scalar with density in [kg/m^3]</p></li></ul><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/CO2Properties/generation.jl#L162-L192" target="_blank" rel="noreferrer">source</a></p>`,11))]),s("details",m,[s("summary",null,[i[21]||(i[21]=s("a",{id:"JutulDarcy.CO2Properties.viscosity_gas_mixture_Davidson1993",href:"#JutulDarcy.CO2Properties.viscosity_gas_mixture_Davidson1993"},[s("span",{class:"jlbinding"},"JutulDarcy.CO2Properties.viscosity_gas_mixture_Davidson1993")],-1)),i[22]||(i[22]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[23]||(i[23]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">viscMixture </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> viscosity_gas_mixture_Davidson1993</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x, M, mu)</span></span></code></pre></div><p>Calculate the viscosity of a gas mixture following Davidson (1993). In principle, valid for any range within which the individual components&#39; viscosities are valid.</p><p>Arguments:</p><ul><li><p>x: Mole fraction of each component</p></li><li><p>M: Molar mass of each component</p></li><li><p>mu: Viscosity of each component in user chosen units</p></li></ul><p>Each input should be a Float64 Vector of length n where n is the total number of components</p><p><strong>Outputs</strong></p><ul><li>viscMixture: Scalar viscosity of the mixture in same units as mu</li></ul><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/CO2Properties/generation.jl#L469-L488" target="_blank" rel="noreferrer">source</a></p>',8))]),i[46]||(i[46]=s("h2",{id:"Relative-permeability-functions",tabindex:"-1"},[t("Relative permeability functions "),s("a",{class:"header-anchor",href:"#Relative-permeability-functions","aria-label":'Permalink to "Relative permeability functions {#Relative-permeability-functions}"'},"​")],-1)),s("details",f,[s("summary",null,[i[24]||(i[24]=s("a",{id:"JutulDarcy.brooks_corey_relperm",href:"#JutulDarcy.brooks_corey_relperm"},[s("span",{class:"jlbinding"},"JutulDarcy.brooks_corey_relperm")],-1)),i[25]||(i[25]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[26]||(i[26]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">brooks_corey_relperm</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(s; n </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 2.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, residual </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, kr_max </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, residual_total </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> residual)</span></span></code></pre></div><p>Evaluate Brooks-Corey relative permeability function at saturation <code>s</code> for exponent <code>n</code> and a given residual and maximum relative permeability value. If considering a two-phase system, the total residual value over both phases should also be passed if the other phase has a non-zero residual value.</p><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/variables/relperm/simple.jl#L176-L183" target="_blank" rel="noreferrer">source</a></p>',3))]),i[47]||(i[47]=s("h2",{id:"CO2-inventory",tabindex:"-1"},[t("CO2 inventory "),s("a",{class:"header-anchor",href:"#CO2-inventory","aria-label":'Permalink to "CO2 inventory {#CO2-inventory}"'},"​")],-1)),s("details",E,[s("summary",null,[i[27]||(i[27]=s("a",{id:"JutulDarcy.co2_inventory",href:"#JutulDarcy.co2_inventory"},[s("span",{class:"jlbinding"},"JutulDarcy.co2_inventory")],-1)),i[28]||(i[28]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[29]||(i[29]=a(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">inventory </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> co2_inventory</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, ws, states, t; cells </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> missing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, co2_name </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;CO2&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">inventory </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> co2_inventory</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, result</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ReservoirSimResult</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; cells </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> missing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Compute CO2 inventory for each step for a given <code>model</code>, well results <code>ws</code> and reporting times t. If provided, the keyword argument <code>cells</code> will compute inventory inside the region defined by the cells, and let any additional CO2 be categorized as &quot;outside region&quot;.</p><p>The inventory will be a Vector of Dicts where each entry contains a breakdown of the status of the CO2 at that time, including residual and dissolution trapping.</p><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/utils.jl#L1851-L1862" target="_blank" rel="noreferrer">source</a></p>`,4))]),s("details",v,[s("summary",null,[i[30]||(i[30]=s("a",{id:"JutulDarcy.plot_co2_inventory",href:"#JutulDarcy.plot_co2_inventory"},[s("span",{class:"jlbinding"},"JutulDarcy.plot_co2_inventory")],-1)),i[31]||(i[31]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[32]||(i[32]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_co2_inventory</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t, inventory, plot_type </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :stack</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Plots the CO2 inventory over time or steps, with options for stacked or line plots. <code>inventory</code> is the output from <code>co2_inventory</code> while <code>t</code> can either be omitted, be a list of reporting time in seconds or a index list of steps where the solution is given.</p><p><strong>Arguments</strong></p><ul><li><p><code>t</code>: A vector representing time or steps. If <code>t</code> is of type <code>Float64</code>, it is assumed to represent time in seconds and will be converted to years.</p></li><li><p><code>inventory</code>: A vector of dictionaries, where each dictionary contains CO2 mass data for different categories (e.g., <code>:dissolved</code>, <code>:mobile</code>, <code>:residual</code>, etc.).</p></li><li><p><code>plot_type</code>: (Optional) A symbol specifying the type of plot. Can be <code>:stack</code> for stacked plots or <code>:lines</code> for line plots. Default is <code>:stack</code>.</p></li></ul><p><strong>Notes</strong></p><p>This function is only available if Makie is loaded (through for example GLMakie or CairoMakie)</p><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/ext.jl#L201-L222" target="_blank" rel="noreferrer">source</a></p>',7))]),i[48]||(i[48]=s("h2",{id:"API-utilities",tabindex:"-1"},[t("API utilities "),s("a",{class:"header-anchor",href:"#API-utilities","aria-label":'Permalink to "API utilities {#API-utilities}"'},"​")],-1)),s("details",C,[s("summary",null,[i[33]||(i[33]=s("a",{id:"JutulDarcy.reservoir_model",href:"#JutulDarcy.reservoir_model"},[s("span",{class:"jlbinding"},"JutulDarcy.reservoir_model")],-1)),i[34]||(i[34]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[35]||(i[35]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">reservoir_model</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model)</span></span></code></pre></div><p>Get the reservoir model from a <code>MultiModel</code> or return the model itself if it is not a <code>MultiModel</code>.</p><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/utils.jl#L1-L6" target="_blank" rel="noreferrer">source</a></p>',3))]),i[49]||(i[49]=s("div",{class:"warning custom-block"},[s("p",{class:"custom-block-title"},"Missing docstring."),s("p",null,[t("Missing docstring for "),s("code",null,"JutulDarcy.reservoir_storage"),t(". Check Documenter's build log for details.")])],-1)),s("details",F,[s("summary",null,[i[36]||(i[36]=s("a",{id:"JutulDarcy.well_symbols",href:"#JutulDarcy.well_symbols"},[s("span",{class:"jlbinding"},"JutulDarcy.well_symbols")],-1)),i[37]||(i[37]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[38]||(i[38]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">well_symbols</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">MultiModel</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Get the keys of a <code>MultiModel</code> models that correspond to well models.</p><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/utils.jl#L1314-L1318" target="_blank" rel="noreferrer">source</a></p>',3))]),i[50]||(i[50]=s("h2",{id:"Adjoints-and-gradients",tabindex:"-1"},[t("Adjoints and gradients "),s("a",{class:"header-anchor",href:"#Adjoints-and-gradients","aria-label":'Permalink to "Adjoints and gradients {#Adjoints-and-gradients}"'},"​")],-1)),s("details",j,[s("summary",null,[i[39]||(i[39]=s("a",{id:"JutulDarcy.reservoir_sensitivities",href:"#JutulDarcy.reservoir_sensitivities"},[s("span",{class:"jlbinding"},"JutulDarcy.reservoir_sensitivities")],-1)),i[40]||(i[40]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[41]||(i[41]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">result, sens </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reservoir_sensitivities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">JutulCase</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, objective</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; sim_arg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> NamedTuple</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), kwarg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Simulate a case and calculate parameter sensitivities with respect to an objective function on the form:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>obj(model, state, dt_n, n, forces_for_step_n)</span></span></code></pre></div><p>The objective is summed up for all steps.</p><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/gradients/gradients.jl#L3-L13" target="_blank" rel="noreferrer">source</a></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">reservoir_sensitivities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">JutulCase</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, rsr</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ReservoirSimResult</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, objective</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; kwarg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Calculate parameter sensitivities with respect to an objective function on the form for a case and a simulation result from that case. The objective function is on the form:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>obj(model, state, dt_n, n, forces_for_step_n)</span></span></code></pre></div><p>The objective is summed up for all steps.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">reservoir_sensitivities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case, rsr, objective; kwarg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/gradients/gradients.jl#L19-L31" target="_blank" rel="noreferrer">source</a></p>',11))]),s("details",_,[s("summary",null,[i[42]||(i[42]=s("a",{id:"JutulDarcy.well_mismatch",href:"#JutulDarcy.well_mismatch"},[s("span",{class:"jlbinding"},"JutulDarcy.well_mismatch")],-1)),i[43]||(i[43]=t()),l(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[44]||(i[44]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">well_mismatch</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(qoi, wells, model_f, states_f, model_c, state_c, dt, step_no, forces; </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">keyword arguments</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Compute well mismatch for a set of qoi&#39;s (well targets) and a set of well symbols.</p><p><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d8fcf98e8e31c6678ff8faddf0059bb49d432184/src/gradients/objectives.jl#L67-L71" target="_blank" rel="noreferrer">source</a></p>',3))])])}const T=n(h,[["render",D]]);export{B as __pageData,T as default};