import{_ as i,o as a,c as e,a5 as n}from"./chunks/framework.K4zi3f8j.js";const c=JSON.parse('{"title":"GPU support","description":"","frontmatter":{},"headers":[],"relativePath":"man/advanced/gpu.md","filePath":"man/advanced/gpu.md","lastUpdated":null}'),t={name:"man/advanced/gpu.md"};function l(h,s,p,r,k,o){return a(),e("div",null,s[0]||(s[0]=[n(`<h1 id="GPU-support" tabindex="-1">GPU support <a class="header-anchor" href="#GPU-support" aria-label="Permalink to &quot;GPU support {#GPU-support}&quot;">​</a></h1><p>JutulDarcy includes experimental support for running linear solves on the GPU. For many simulations, the linear systems are the most compute-intensive part and a natural choice for acceleration. At the moment, the support is limited to CUDA GPUs through <a href="https://github.com/JuliaGPU/CUDA.jl" target="_blank" rel="noreferrer">CUDA.jl</a>. For the most efficient CPR preconditioner, <a href="https://github.com/JuliaGPU/AMGX.jl" target="_blank" rel="noreferrer">AMGX.jl</a> is required which is currently limited to Linux systems. Windows users may have luck by running Julia inside <a href="https://learn.microsoft.com/en-us/windows/wsl/install" target="_blank" rel="noreferrer">WSL</a>.</p><h2 id="How-to-use" tabindex="-1">How to use <a class="header-anchor" href="#How-to-use" aria-label="Permalink to &quot;How to use {#How-to-use}&quot;">​</a></h2><p>If you have installed JutulDarcy, you should start by adding the CUDA and optionally the AMGX packages using the package manager:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;CUDA&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Requires a CUDA-capable GPU</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;AMGX&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Requires CUDA + Linux</span></span></code></pre></div><p>Once the packages have been added to the same environment as JutulDarcy, you can load them to enable GPU support. Let us grab the first ten steps of the EGG benchmark model:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Jutul, JutulDarcy</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dpth </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">GeoEnergyIO</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">test_input_file_path</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;EGG&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;EGG.DATA&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">case </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_case_from_data_file</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(dpth)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">case </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> case[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span></code></pre></div><h3 id="Running-on-CPU" tabindex="-1">Running on CPU <a class="header-anchor" href="#Running-on-CPU" aria-label="Permalink to &quot;Running on CPU {#Running-on-CPU}&quot;">​</a></h3><p>If we wanted to run this on CPU we would simply call <code>simulate_reservoir</code>:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">result_cpu </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case);</span></span></code></pre></div><h3 id="Running-on-GPU-with-block-ILU(0)" tabindex="-1">Running on GPU with block ILU(0) <a class="header-anchor" href="#Running-on-GPU-with-block-ILU(0)" aria-label="Permalink to &quot;Running on GPU with block ILU(0) {#Running-on-GPU-with-block-ILU(0)}&quot;">​</a></h3><p>If we now load <code>CUDA</code> we can run the same simulation using the CUDA-accelerated linear solver. By itself, CUDA only supports the ILU(0) preconditioner. JutulDarcy will automatically pick this preconditioner when CUDA is requested without AMGX, but we write it explicitly here:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CUDA</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">result_ilu0_cuda </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case, linear_solver_backend </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :cuda</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, precond </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :ilu0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">);</span></span></code></pre></div><h3 id="Running-on-GPU-with-CPR-AMGX-ILU(0)" tabindex="-1">Running on GPU with CPR AMGX-ILU(0) <a class="header-anchor" href="#Running-on-GPU-with-CPR-AMGX-ILU(0)" aria-label="Permalink to &quot;Running on GPU with CPR AMGX-ILU(0) {#Running-on-GPU-with-CPR-AMGX-ILU(0)}&quot;">​</a></h3><p>Loading the AMGX package makes a pure GPU-based two-stage CPR available. Again, we are explicit in requesting CPR, but if both <code>CUDA</code> and <code>AMGX</code> are available and functional this is redundant:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> AMGX</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">result_amgx_cuda </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case, linear_solver_backend </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :cuda</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, precond </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :cpr</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">);</span></span></code></pre></div><p>In short, load <code>AMGX</code> and <code>CUDA</code> and run <code>simulate_reservoir(case, linear_solver_backend = :cuda)</code> to get GPU results. The EGG model is quite small, so if you want to see significant performance increases, a larger case will be necessary. <code>AMGX</code> also contains a large number of options that can be configured for advanced users.</p><h2 id="Technical-details-and-limitations" tabindex="-1">Technical details and limitations <a class="header-anchor" href="#Technical-details-and-limitations" aria-label="Permalink to &quot;Technical details and limitations {#Technical-details-and-limitations}&quot;">​</a></h2><p>The GPU implementation relies on assembly on CPU and pinned memory to transfer onto the CPU. This means that the performance can be significantly improved by launching Julia with multiple threads to speed up the non-GPU parts of the code. AMGX is currently single-GPU only and does not work with MPI. To make use of lower precision, specify <code>Float32</code> in the <code>float_type</code> argument to the linear solver. Additional arguments to <code>AMGX</code> can also be specified this way. For example, we can solve using aggregation AMG in single precision by doing the following:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">simulate_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    linear_solver_backend </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :cuda</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    linear_solver_arg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        float_type </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Float32,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        algorithm </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;AGGREGATION&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        selector </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;SIZE_8&quot;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        )</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    )</span></span></code></pre></div><div class="warning custom-block"><p class="custom-block-title">Experimental status</p><p>Multiple successive runs with different <code>AMGX</code> instances have resulted in crashes when old instances are garbage collected. This part of the code is still considered experimental, with contributions welcome if you are using it.</p></div>`,21)]))}const u=i(t,[["render",l]]);export{c as __pageData,u as default};