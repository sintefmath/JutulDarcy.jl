import{_ as t,c as n,j as a,a as e,a5 as i,o as p}from"./chunks/framework.B4o_rkVA.js";const l="/JutulDarcy.jl/v0.2.33/assets/eavkrki.BIVOWHXu.jpeg",k="/JutulDarcy.jl/v0.2.33/assets/cmgjphh.COf0RfYU.jpeg",c=JSON.parse('{"title":"Validation of equation-of-state compositiona simulator","description":"","frontmatter":{},"headers":[],"relativePath":"examples/validation_compositional.md","filePath":"examples/validation_compositional.md","lastUpdated":null}'),y={name:"examples/validation_compositional.md"},r={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},h={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"0.988ex",height:"1.065ex",role:"img",focusable:"false",viewBox:"0 -320.9 436.6 470.9","aria-hidden":"true"},o={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},w={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"0.988ex",height:"1.065ex",role:"img",focusable:"false",viewBox:"0 -320.9 436.6 470.9","aria-hidden":"true"};function u(d,s,m,g,b,S){return p(),n("div",null,[s[9]||(s[9]=a("h1",{id:"Validation-of-equation-of-state-compositiona-simulator",tabindex:"-1"},[e("Validation of equation-of-state compositiona simulator "),a("a",{class:"header-anchor",href:"#Validation-of-equation-of-state-compositiona-simulator","aria-label":'Permalink to "Validation of equation-of-state compositiona simulator {#Validation-of-equation-of-state-compositiona-simulator}"'},"​")],-1)),s[10]||(s[10]=a("p",null,"This example solves a 1D two-phase, three component miscible displacement problem and compares against existing simulators (E300, AD-GPRS) to verify correctness.",-1)),s[11]||(s[11]=a("p",null,"The case is loaded from an input file that can be run in other simulators. For convenience, we provide solutions from the other simulators as a binary file to perform a comparison without having to run and convert results from other the simulators.",-1)),a("p",null,[s[4]||(s[4]=e("This case is a small compositional problem inspired by the examples in Voskov et al (JPSE, 2012). A 1D reservoir of 1,000 meters length is discretized into 1,000 cells. The model initially contains a mixture made up of 0.6 parts C10, 0.1 parts CO2, and 0.3 parts C1 by moles at 150 degrees C and 75 bar pressure. Wells are placed in the leftmost and rightmost cells of the domain, with the leftmost well injecting pure CO")),a("mjx-container",r,[(p(),n("svg",h,s[0]||(s[0]=[i('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="msub"><g data-mml-node="mi"></g><g data-mml-node="mn" transform="translate(33,-150) scale(0.707)"><path data-c="32" d="M109 429Q82 429 66 447T50 491Q50 562 103 614T235 666Q326 666 387 610T449 465Q449 422 429 383T381 315T301 241Q265 210 201 149L142 93L218 92Q375 92 385 97Q392 99 409 186V189H449V186Q448 183 436 95T421 3V0H50V19V31Q50 38 56 46T86 81Q115 113 136 137Q145 147 170 174T204 211T233 244T261 278T284 308T305 340T320 369T333 401T340 431T343 464Q343 527 309 573T212 619Q179 619 154 602T119 569T109 550Q109 549 114 549Q132 549 151 535T170 489Q170 464 154 447T109 429Z" style="stroke-width:3;"></path></g></g></g></g>',1)]))),s[1]||(s[1]=a("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[a("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[a("msub",null,[a("mi"),a("mn",null,"2")])])],-1))]),s[5]||(s[5]=e(" at a fixed bottom-hole pressure of 100 bar and the other well producing at 50 bar. The model is isothermal and contains a phase transition from the initial two-phase mixture to single-phase gas as injected CO")),a("mjx-container",o,[(p(),n("svg",w,s[2]||(s[2]=[i('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="msub"><g data-mml-node="mi"></g><g data-mml-node="mn" transform="translate(33,-150) scale(0.707)"><path data-c="32" d="M109 429Q82 429 66 447T50 491Q50 562 103 614T235 666Q326 666 387 610T449 465Q449 422 429 383T381 315T301 241Q265 210 201 149L142 93L218 92Q375 92 385 97Q392 99 409 186V189H449V186Q448 183 436 95T421 3V0H50V19V31Q50 38 56 46T86 81Q115 113 136 137Q145 147 170 174T204 211T233 244T261 278T284 308T305 340T320 369T333 401T340 431T343 464Q343 527 309 573T212 619Q179 619 154 602T119 569T109 550Q109 549 114 549Q132 549 151 535T170 489Q170 464 154 447T109 429Z" style="stroke-width:3;"></path></g></g></g></g>',1)]))),s[3]||(s[3]=a("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[a("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[a("msub",null,[a("mi"),a("mn",null,"2")])])],-1))]),s[6]||(s[6]=e(" eventually displaces the resident fluids. For further details on this setup, see Møyner and Tchelepi (SPE J. 2018) [")),s[7]||(s[7]=a("a",{href:"/JutulDarcy.jl/v0.2.33/extras/refs#moyner_tchelepi_2018"},"8",-1)),s[8]||(s[8]=e("]."))]),s[12]||(s[12]=i(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Jutul</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GLMakie</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dpth </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">GeoEnergyIO</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">test_input_file_path</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;SIMPLE_COMP&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">data_path </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> joinpath</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(dpth, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;SIMPLE_COMP.DATA&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">case </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> setup_case_from_data_file</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(data_path)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">result </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> simulate_reservoir</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case, info_level </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ws, states </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> result;</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>\x1B[92;1mJutul:\x1B[0m Simulating 19 years, 32.21 weeks as 500 report steps</span></span>
<span class="line"><span>\x1B[34;1mStep   1/500:\x1B[0m Solving start to 14 minutes, 24 seconds, Δt = 14 minutes, 24 seconds</span></span>
<span class="line"><span>\x1B[34;1mStep   2/500:\x1B[0m Solving 14 minutes, 24 seconds to 43 minutes, 12 seconds, Δt = 28 minutes, 48 seconds</span></span>
<span class="line"><span>\x1B[34;1mStep   3/500:\x1B[0m Solving 43 minutes, 12 seconds to 1 hour, 40.8 minutes, Δt = 57 minutes, 36 seconds</span></span>
<span class="line"><span>\x1B[34;1mStep   4/500:\x1B[0m Solving 1 hour, 40.8 minutes to 2 hours, 38.4 minutes, Δt = 57 minutes, 36 seconds</span></span>
<span class="line"><span>\x1B[34;1mStep   5/500:\x1B[0m Solving 2 hours, 38.4 minutes to 4 hours, 33.6 minutes, Δt = 1 hour, 55.2 minutes</span></span>
<span class="line"><span>\x1B[34;1mStep   6/500:\x1B[0m Solving 4 hours, 33.6 minutes to 8 hours, 24 minutes, Δt = 3 hours, 50.4 minutes</span></span>
<span class="line"><span>\x1B[34;1mStep   7/500:\x1B[0m Solving 8 hours, 24 minutes to 16 hours, 4.8 minutes, Δt = 7 hours, 40.8 minutes</span></span>
<span class="line"><span>\x1B[34;1mStep   8/500:\x1B[0m Solving 16 hours, 4.8 minutes to 1 day, 2.64 hours, Δt = 10 hours, 33.6 minutes</span></span>
<span class="line"><span>\x1B[34;1mStep   9/500:\x1B[0m Solving 1 day, 2.64 hours to 1 day, 23.76 hours, Δt = 21 hours, 7.2 minutes</span></span>
<span class="line"><span>\x1B[34;1mStep  10/500:\x1B[0m Solving 1 day, 23.76 hours to 2 days, 2.64 hours, Δt = 2 hours, 52.8 minutes</span></span>
<span class="line"><span>\x1B[34;1mStep  11/500:\x1B[0m Solving 2 days, 2.64 hours to 3 days, 2.64 hours, Δt = 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  12/500:\x1B[0m Solving 3 days, 2.64 hours to 4 days, 2.64 hours, Δt = 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  13/500:\x1B[0m Solving 4 days, 2.64 hours to 5 days, 2.64 hours, Δt = 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  14/500:\x1B[0m Solving 5 days, 2.64 hours to 1 week, 2.64 hours, Δt = 2 days</span></span>
<span class="line"><span>\x1B[34;1mStep  15/500:\x1B[0m Solving 1 week, 2.64 hours to 1 week, 3.11 days, Δt = 3 days</span></span>
<span class="line"><span>\x1B[34;1mStep  16/500:\x1B[0m Solving 1 week, 3.11 days to 2 weeks, 1.11 day, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  17/500:\x1B[0m Solving 2 weeks, 1.11 day to 2 weeks, 6.11 days, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  18/500:\x1B[0m Solving 2 weeks, 6.11 days to 3 weeks, 4.11 days, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  19/500:\x1B[0m Solving 3 weeks, 4.11 days to 4 weeks, 2.11 days, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  20/500:\x1B[0m Solving 4 weeks, 2.11 days to 5 weeks, 2.64 hours, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  21/500:\x1B[0m Solving 5 weeks, 2.64 hours to 5 weeks, 5.11 days, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  22/500:\x1B[0m Solving 5 weeks, 5.11 days to 6 weeks, 3.11 days, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  23/500:\x1B[0m Solving 6 weeks, 3.11 days to 7 weeks, 1.11 day, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  24/500:\x1B[0m Solving 7 weeks, 1.11 day to 7 weeks, 6.11 days, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  25/500:\x1B[0m Solving 7 weeks, 6.11 days to 8 weeks, 4.11 days, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  26/500:\x1B[0m Solving 8 weeks, 4.11 days to 9 weeks, 2.11 days, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  27/500:\x1B[0m Solving 9 weeks, 2.11 days to 10 weeks, 2.64 hours, Δt = 5 days</span></span>
<span class="line"><span>\x1B[34;1mStep  28/500:\x1B[0m Solving 10 weeks, 2.64 hours to 12 weeks, 1.11 day, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  29/500:\x1B[0m Solving 12 weeks, 1.11 day to 14 weeks, 2.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  30/500:\x1B[0m Solving 14 weeks, 2.11 days to 16 weeks, 3.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  31/500:\x1B[0m Solving 16 weeks, 3.11 days to 18 weeks, 4.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  32/500:\x1B[0m Solving 18 weeks, 4.11 days to 20 weeks, 5.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  33/500:\x1B[0m Solving 20 weeks, 5.11 days to 22 weeks, 6.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  34/500:\x1B[0m Solving 22 weeks, 6.11 days to 25 weeks, 2.64 hours, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  35/500:\x1B[0m Solving 25 weeks, 2.64 hours to 27 weeks, 1.11 day, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  36/500:\x1B[0m Solving 27 weeks, 1.11 day to 29 weeks, 2.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  37/500:\x1B[0m Solving 29 weeks, 2.11 days to 31 weeks, 3.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  38/500:\x1B[0m Solving 31 weeks, 3.11 days to 33 weeks, 4.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  39/500:\x1B[0m Solving 33 weeks, 4.11 days to 35 weeks, 5.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  40/500:\x1B[0m Solving 35 weeks, 5.11 days to 37 weeks, 6.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  41/500:\x1B[0m Solving 37 weeks, 6.11 days to 40 weeks, 2.64 hours, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  42/500:\x1B[0m Solving 40 weeks, 2.64 hours to 42 weeks, 1.11 day, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  43/500:\x1B[0m Solving 42 weeks, 1.11 day to 44 weeks, 2.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  44/500:\x1B[0m Solving 44 weeks, 2.11 days to 46 weeks, 3.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  45/500:\x1B[0m Solving 46 weeks, 3.11 days to 48 weeks, 4.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  46/500:\x1B[0m Solving 48 weeks, 4.11 days to 50 weeks, 5.11 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  47/500:\x1B[0m Solving 50 weeks, 5.11 days to 1 year, 4.867 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  48/500:\x1B[0m Solving 1 year, 4.867 days to 1 year, 2.838 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  49/500:\x1B[0m Solving 1 year, 2.838 weeks to 1 year, 4.981 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  50/500:\x1B[0m Solving 1 year, 4.981 weeks to 1 year, 7.124 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  51/500:\x1B[0m Solving 1 year, 7.124 weeks to 1 year, 9.267 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  52/500:\x1B[0m Solving 1 year, 9.267 weeks to 1 year, 11.41 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  53/500:\x1B[0m Solving 1 year, 11.41 weeks to 1 year, 13.55 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  54/500:\x1B[0m Solving 1 year, 13.55 weeks to 1 year, 15.7 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  55/500:\x1B[0m Solving 1 year, 15.7 weeks to 1 year, 17.84 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  56/500:\x1B[0m Solving 1 year, 17.84 weeks to 1 year, 19.98 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  57/500:\x1B[0m Solving 1 year, 19.98 weeks to 1 year, 22.12 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  58/500:\x1B[0m Solving 1 year, 22.12 weeks to 1 year, 24.27 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  59/500:\x1B[0m Solving 1 year, 24.27 weeks to 1 year, 26.41 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  60/500:\x1B[0m Solving 1 year, 26.41 weeks to 1 year, 28.55 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  61/500:\x1B[0m Solving 1 year, 28.55 weeks to 1 year, 30.7 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  62/500:\x1B[0m Solving 1 year, 30.7 weeks to 1 year, 32.84 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  63/500:\x1B[0m Solving 1 year, 32.84 weeks to 1 year, 34.98 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  64/500:\x1B[0m Solving 1 year, 34.98 weeks to 1 year, 37.12 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  65/500:\x1B[0m Solving 1 year, 37.12 weeks to 1 year, 39.27 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  66/500:\x1B[0m Solving 1 year, 39.27 weeks to 1 year, 41.41 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  67/500:\x1B[0m Solving 1 year, 41.41 weeks to 1 year, 43.55 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  68/500:\x1B[0m Solving 1 year, 43.55 weeks to 1 year, 45.7 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  69/500:\x1B[0m Solving 1 year, 45.7 weeks to 1 year, 47.84 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  70/500:\x1B[0m Solving 1 year, 47.84 weeks to 1 year, 49.98 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  71/500:\x1B[0m Solving 1 year, 49.98 weeks to 1 year, 52.12 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  72/500:\x1B[0m Solving 1 year, 52.12 weeks to 2 years, 2.089 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  73/500:\x1B[0m Solving 2 years, 2.089 weeks to 2 years, 4.232 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  74/500:\x1B[0m Solving 2 years, 4.232 weeks to 2 years, 6.375 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  75/500:\x1B[0m Solving 2 years, 6.375 weeks to 2 years, 8.518 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  76/500:\x1B[0m Solving 2 years, 8.518 weeks to 2 years, 10.66 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  77/500:\x1B[0m Solving 2 years, 10.66 weeks to 2 years, 12.8 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  78/500:\x1B[0m Solving 2 years, 12.8 weeks to 2 years, 14.95 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  79/500:\x1B[0m Solving 2 years, 14.95 weeks to 2 years, 17.09 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  80/500:\x1B[0m Solving 2 years, 17.09 weeks to 2 years, 19.23 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  81/500:\x1B[0m Solving 2 years, 19.23 weeks to 2 years, 21.38 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  82/500:\x1B[0m Solving 2 years, 21.38 weeks to 2 years, 23.52 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  83/500:\x1B[0m Solving 2 years, 23.52 weeks to 2 years, 25.66 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  84/500:\x1B[0m Solving 2 years, 25.66 weeks to 2 years, 27.8 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  85/500:\x1B[0m Solving 2 years, 27.8 weeks to 2 years, 29.95 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  86/500:\x1B[0m Solving 2 years, 29.95 weeks to 2 years, 32.09 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  87/500:\x1B[0m Solving 2 years, 32.09 weeks to 2 years, 34.23 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  88/500:\x1B[0m Solving 2 years, 34.23 weeks to 2 years, 36.38 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  89/500:\x1B[0m Solving 2 years, 36.38 weeks to 2 years, 38.52 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  90/500:\x1B[0m Solving 2 years, 38.52 weeks to 2 years, 40.66 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  91/500:\x1B[0m Solving 2 years, 40.66 weeks to 2 years, 42.8 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  92/500:\x1B[0m Solving 2 years, 42.8 weeks to 2 years, 44.95 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  93/500:\x1B[0m Solving 2 years, 44.95 weeks to 2 years, 47.09 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  94/500:\x1B[0m Solving 2 years, 47.09 weeks to 2 years, 49.23 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  95/500:\x1B[0m Solving 2 years, 49.23 weeks to 2 years, 51.38 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  96/500:\x1B[0m Solving 2 years, 51.38 weeks to 3 years, 1.34 week, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  97/500:\x1B[0m Solving 3 years, 1.34 week to 3 years, 3.483 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  98/500:\x1B[0m Solving 3 years, 3.483 weeks to 3 years, 5.626 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep  99/500:\x1B[0m Solving 3 years, 5.626 weeks to 3 years, 7.769 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 100/500:\x1B[0m Solving 3 years, 7.769 weeks to 3 years, 9.912 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 101/500:\x1B[0m Solving 3 years, 9.912 weeks to 3 years, 12.05 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 102/500:\x1B[0m Solving 3 years, 12.05 weeks to 3 years, 14.2 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 103/500:\x1B[0m Solving 3 years, 14.2 weeks to 3 years, 16.34 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 104/500:\x1B[0m Solving 3 years, 16.34 weeks to 3 years, 18.48 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 105/500:\x1B[0m Solving 3 years, 18.48 weeks to 3 years, 20.63 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 106/500:\x1B[0m Solving 3 years, 20.63 weeks to 3 years, 22.77 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 107/500:\x1B[0m Solving 3 years, 22.77 weeks to 3 years, 24.91 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 108/500:\x1B[0m Solving 3 years, 24.91 weeks to 3 years, 27.05 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 109/500:\x1B[0m Solving 3 years, 27.05 weeks to 3 years, 29.2 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 110/500:\x1B[0m Solving 3 years, 29.2 weeks to 3 years, 31.34 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 111/500:\x1B[0m Solving 3 years, 31.34 weeks to 3 years, 33.48 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 112/500:\x1B[0m Solving 3 years, 33.48 weeks to 3 years, 35.63 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 113/500:\x1B[0m Solving 3 years, 35.63 weeks to 3 years, 37.77 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 114/500:\x1B[0m Solving 3 years, 37.77 weeks to 3 years, 39.91 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 115/500:\x1B[0m Solving 3 years, 39.91 weeks to 3 years, 42.05 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 116/500:\x1B[0m Solving 3 years, 42.05 weeks to 3 years, 44.2 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 117/500:\x1B[0m Solving 3 years, 44.2 weeks to 3 years, 46.34 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 118/500:\x1B[0m Solving 3 years, 46.34 weeks to 3 years, 48.48 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 119/500:\x1B[0m Solving 3 years, 48.48 weeks to 3 years, 50.63 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 120/500:\x1B[0m Solving 3 years, 50.63 weeks to 4 years, 4.14 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 121/500:\x1B[0m Solving 4 years, 4.14 days to 4 years, 2.734 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 122/500:\x1B[0m Solving 4 years, 2.734 weeks to 4 years, 4.877 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 123/500:\x1B[0m Solving 4 years, 4.877 weeks to 4 years, 7.02 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 124/500:\x1B[0m Solving 4 years, 7.02 weeks to 4 years, 9.163 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 125/500:\x1B[0m Solving 4 years, 9.163 weeks to 4 years, 11.31 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 126/500:\x1B[0m Solving 4 years, 11.31 weeks to 4 years, 13.45 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 127/500:\x1B[0m Solving 4 years, 13.45 weeks to 4 years, 15.59 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 128/500:\x1B[0m Solving 4 years, 15.59 weeks to 4 years, 17.73 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 129/500:\x1B[0m Solving 4 years, 17.73 weeks to 4 years, 19.88 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 130/500:\x1B[0m Solving 4 years, 19.88 weeks to 4 years, 22.02 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 131/500:\x1B[0m Solving 4 years, 22.02 weeks to 4 years, 24.16 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 132/500:\x1B[0m Solving 4 years, 24.16 weeks to 4 years, 26.31 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 133/500:\x1B[0m Solving 4 years, 26.31 weeks to 4 years, 28.45 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 134/500:\x1B[0m Solving 4 years, 28.45 weeks to 4 years, 30.59 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 135/500:\x1B[0m Solving 4 years, 30.59 weeks to 4 years, 32.73 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 136/500:\x1B[0m Solving 4 years, 32.73 weeks to 4 years, 34.88 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 137/500:\x1B[0m Solving 4 years, 34.88 weeks to 4 years, 37.02 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 138/500:\x1B[0m Solving 4 years, 37.02 weeks to 4 years, 39.16 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 139/500:\x1B[0m Solving 4 years, 39.16 weeks to 4 years, 41.31 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 140/500:\x1B[0m Solving 4 years, 41.31 weeks to 4 years, 43.45 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 141/500:\x1B[0m Solving 4 years, 43.45 weeks to 4 years, 45.59 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 142/500:\x1B[0m Solving 4 years, 45.59 weeks to 4 years, 47.73 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 143/500:\x1B[0m Solving 4 years, 47.73 weeks to 4 years, 49.88 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 144/500:\x1B[0m Solving 4 years, 49.88 weeks to 4 years, 52.02 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 145/500:\x1B[0m Solving 4 years, 52.02 weeks to 5 years, 1.985 week, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 146/500:\x1B[0m Solving 5 years, 1.985 week to 5 years, 4.128 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 147/500:\x1B[0m Solving 5 years, 4.128 weeks to 5 years, 6.271 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 148/500:\x1B[0m Solving 5 years, 6.271 weeks to 5 years, 8.414 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 149/500:\x1B[0m Solving 5 years, 8.414 weeks to 5 years, 10.56 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 150/500:\x1B[0m Solving 5 years, 10.56 weeks to 5 years, 12.7 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 151/500:\x1B[0m Solving 5 years, 12.7 weeks to 5 years, 14.84 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 152/500:\x1B[0m Solving 5 years, 14.84 weeks to 5 years, 16.99 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 153/500:\x1B[0m Solving 5 years, 16.99 weeks to 5 years, 19.13 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 154/500:\x1B[0m Solving 5 years, 19.13 weeks to 5 years, 21.27 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 155/500:\x1B[0m Solving 5 years, 21.27 weeks to 5 years, 23.41 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 156/500:\x1B[0m Solving 5 years, 23.41 weeks to 5 years, 25.56 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 157/500:\x1B[0m Solving 5 years, 25.56 weeks to 5 years, 27.7 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 158/500:\x1B[0m Solving 5 years, 27.7 weeks to 5 years, 29.84 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 159/500:\x1B[0m Solving 5 years, 29.84 weeks to 5 years, 31.99 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 160/500:\x1B[0m Solving 5 years, 31.99 weeks to 5 years, 34.13 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 161/500:\x1B[0m Solving 5 years, 34.13 weeks to 5 years, 36.27 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 162/500:\x1B[0m Solving 5 years, 36.27 weeks to 5 years, 38.41 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 163/500:\x1B[0m Solving 5 years, 38.41 weeks to 5 years, 40.56 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 164/500:\x1B[0m Solving 5 years, 40.56 weeks to 5 years, 42.7 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 165/500:\x1B[0m Solving 5 years, 42.7 weeks to 5 years, 44.84 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 166/500:\x1B[0m Solving 5 years, 44.84 weeks to 5 years, 46.99 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 167/500:\x1B[0m Solving 5 years, 46.99 weeks to 5 years, 49.13 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 168/500:\x1B[0m Solving 5 years, 49.13 weeks to 5 years, 51.27 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 169/500:\x1B[0m Solving 5 years, 51.27 weeks to 6 years, 1.236 week, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 170/500:\x1B[0m Solving 6 years, 1.236 week to 6 years, 3.379 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 171/500:\x1B[0m Solving 6 years, 3.379 weeks to 6 years, 5.522 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 172/500:\x1B[0m Solving 6 years, 5.522 weeks to 6 years, 7.665 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 173/500:\x1B[0m Solving 6 years, 7.665 weeks to 6 years, 9.808 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 174/500:\x1B[0m Solving 6 years, 9.808 weeks to 6 years, 11.95 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 175/500:\x1B[0m Solving 6 years, 11.95 weeks to 6 years, 14.09 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 176/500:\x1B[0m Solving 6 years, 14.09 weeks to 6 years, 16.24 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 177/500:\x1B[0m Solving 6 years, 16.24 weeks to 6 years, 18.38 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 178/500:\x1B[0m Solving 6 years, 18.38 weeks to 6 years, 20.52 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 179/500:\x1B[0m Solving 6 years, 20.52 weeks to 6 years, 22.66 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 180/500:\x1B[0m Solving 6 years, 22.66 weeks to 6 years, 24.81 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 181/500:\x1B[0m Solving 6 years, 24.81 weeks to 6 years, 26.95 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 182/500:\x1B[0m Solving 6 years, 26.95 weeks to 6 years, 29.09 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 183/500:\x1B[0m Solving 6 years, 29.09 weeks to 6 years, 31.24 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 184/500:\x1B[0m Solving 6 years, 31.24 weeks to 6 years, 33.38 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 185/500:\x1B[0m Solving 6 years, 33.38 weeks to 6 years, 35.52 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 186/500:\x1B[0m Solving 6 years, 35.52 weeks to 6 years, 37.66 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 187/500:\x1B[0m Solving 6 years, 37.66 weeks to 6 years, 39.81 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 188/500:\x1B[0m Solving 6 years, 39.81 weeks to 6 years, 41.95 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 189/500:\x1B[0m Solving 6 years, 41.95 weeks to 6 years, 44.09 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 190/500:\x1B[0m Solving 6 years, 44.09 weeks to 6 years, 46.24 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 191/500:\x1B[0m Solving 6 years, 46.24 weeks to 6 years, 48.38 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 192/500:\x1B[0m Solving 6 years, 48.38 weeks to 6 years, 50.52 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 193/500:\x1B[0m Solving 6 years, 50.52 weeks to 7 years, 3.413 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 194/500:\x1B[0m Solving 7 years, 3.413 days to 7 years, 2.63 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 195/500:\x1B[0m Solving 7 years, 2.63 weeks to 7 years, 4.773 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 196/500:\x1B[0m Solving 7 years, 4.773 weeks to 7 years, 6.916 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 197/500:\x1B[0m Solving 7 years, 6.916 weeks to 7 years, 9.059 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 198/500:\x1B[0m Solving 7 years, 9.059 weeks to 7 years, 11.2 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 199/500:\x1B[0m Solving 7 years, 11.2 weeks to 7 years, 13.34 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 200/500:\x1B[0m Solving 7 years, 13.34 weeks to 7 years, 15.49 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 201/500:\x1B[0m Solving 7 years, 15.49 weeks to 7 years, 17.63 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 202/500:\x1B[0m Solving 7 years, 17.63 weeks to 7 years, 19.77 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 203/500:\x1B[0m Solving 7 years, 19.77 weeks to 7 years, 21.92 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 204/500:\x1B[0m Solving 7 years, 21.92 weeks to 7 years, 24.06 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 205/500:\x1B[0m Solving 7 years, 24.06 weeks to 7 years, 26.2 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 206/500:\x1B[0m Solving 7 years, 26.2 weeks to 7 years, 28.34 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 207/500:\x1B[0m Solving 7 years, 28.34 weeks to 7 years, 30.49 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 208/500:\x1B[0m Solving 7 years, 30.49 weeks to 7 years, 32.63 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 209/500:\x1B[0m Solving 7 years, 32.63 weeks to 7 years, 34.77 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 210/500:\x1B[0m Solving 7 years, 34.77 weeks to 7 years, 36.92 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 211/500:\x1B[0m Solving 7 years, 36.92 weeks to 7 years, 39.06 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 212/500:\x1B[0m Solving 7 years, 39.06 weeks to 7 years, 41.2 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 213/500:\x1B[0m Solving 7 years, 41.2 weeks to 7 years, 43.34 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 214/500:\x1B[0m Solving 7 years, 43.34 weeks to 7 years, 45.49 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 215/500:\x1B[0m Solving 7 years, 45.49 weeks to 7 years, 47.63 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 216/500:\x1B[0m Solving 7 years, 47.63 weeks to 7 years, 49.77 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 217/500:\x1B[0m Solving 7 years, 49.77 weeks to 7 years, 51.92 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 218/500:\x1B[0m Solving 7 years, 51.92 weeks to 8 years, 1.881 week, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 219/500:\x1B[0m Solving 8 years, 1.881 week to 8 years, 4.024 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 220/500:\x1B[0m Solving 8 years, 4.024 weeks to 8 years, 6.167 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 221/500:\x1B[0m Solving 8 years, 6.167 weeks to 8 years, 8.31 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 222/500:\x1B[0m Solving 8 years, 8.31 weeks to 8 years, 10.45 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 223/500:\x1B[0m Solving 8 years, 10.45 weeks to 8 years, 12.6 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 224/500:\x1B[0m Solving 8 years, 12.6 weeks to 8 years, 14.74 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 225/500:\x1B[0m Solving 8 years, 14.74 weeks to 8 years, 16.88 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 226/500:\x1B[0m Solving 8 years, 16.88 weeks to 8 years, 19.02 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 227/500:\x1B[0m Solving 8 years, 19.02 weeks to 8 years, 21.17 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 228/500:\x1B[0m Solving 8 years, 21.17 weeks to 8 years, 23.31 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 229/500:\x1B[0m Solving 8 years, 23.31 weeks to 8 years, 25.45 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 230/500:\x1B[0m Solving 8 years, 25.45 weeks to 8 years, 27.6 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 231/500:\x1B[0m Solving 8 years, 27.6 weeks to 8 years, 29.74 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 232/500:\x1B[0m Solving 8 years, 29.74 weeks to 8 years, 31.88 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 233/500:\x1B[0m Solving 8 years, 31.88 weeks to 8 years, 34.02 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 234/500:\x1B[0m Solving 8 years, 34.02 weeks to 8 years, 36.17 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 235/500:\x1B[0m Solving 8 years, 36.17 weeks to 8 years, 38.31 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 236/500:\x1B[0m Solving 8 years, 38.31 weeks to 8 years, 40.45 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 237/500:\x1B[0m Solving 8 years, 40.45 weeks to 8 years, 42.6 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 238/500:\x1B[0m Solving 8 years, 42.6 weeks to 8 years, 44.74 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 239/500:\x1B[0m Solving 8 years, 44.74 weeks to 8 years, 46.88 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 240/500:\x1B[0m Solving 8 years, 46.88 weeks to 8 years, 49.02 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 241/500:\x1B[0m Solving 8 years, 49.02 weeks to 8 years, 51.17 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 242/500:\x1B[0m Solving 8 years, 51.17 weeks to 9 years, 1.133 week, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 243/500:\x1B[0m Solving 9 years, 1.133 week to 9 years, 3.275 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 244/500:\x1B[0m Solving 9 years, 3.275 weeks to 9 years, 5.418 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 245/500:\x1B[0m Solving 9 years, 5.418 weeks to 9 years, 7.561 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 246/500:\x1B[0m Solving 9 years, 7.561 weeks to 9 years, 9.704 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 247/500:\x1B[0m Solving 9 years, 9.704 weeks to 9 years, 11.85 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 248/500:\x1B[0m Solving 9 years, 11.85 weeks to 9 years, 13.99 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 249/500:\x1B[0m Solving 9 years, 13.99 weeks to 9 years, 16.13 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 250/500:\x1B[0m Solving 9 years, 16.13 weeks to 9 years, 18.28 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 251/500:\x1B[0m Solving 9 years, 18.28 weeks to 9 years, 20.42 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 252/500:\x1B[0m Solving 9 years, 20.42 weeks to 9 years, 22.56 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 253/500:\x1B[0m Solving 9 years, 22.56 weeks to 9 years, 24.7 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 254/500:\x1B[0m Solving 9 years, 24.7 weeks to 9 years, 26.85 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 255/500:\x1B[0m Solving 9 years, 26.85 weeks to 9 years, 28.99 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 256/500:\x1B[0m Solving 9 years, 28.99 weeks to 9 years, 31.13 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 257/500:\x1B[0m Solving 9 years, 31.13 weeks to 9 years, 33.28 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 258/500:\x1B[0m Solving 9 years, 33.28 weeks to 9 years, 35.42 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 259/500:\x1B[0m Solving 9 years, 35.42 weeks to 9 years, 37.56 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 260/500:\x1B[0m Solving 9 years, 37.56 weeks to 9 years, 39.7 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 261/500:\x1B[0m Solving 9 years, 39.7 weeks to 9 years, 41.85 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 262/500:\x1B[0m Solving 9 years, 41.85 weeks to 9 years, 43.99 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 263/500:\x1B[0m Solving 9 years, 43.99 weeks to 9 years, 46.13 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 264/500:\x1B[0m Solving 9 years, 46.13 weeks to 9 years, 48.28 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 265/500:\x1B[0m Solving 9 years, 48.28 weeks to 9 years, 50.42 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 266/500:\x1B[0m Solving 9 years, 50.42 weeks to 10 years, 2.685 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 267/500:\x1B[0m Solving 10 years, 2.685 days to 10 years, 2.526 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 268/500:\x1B[0m Solving 10 years, 2.526 weeks to 10 years, 4.669 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 269/500:\x1B[0m Solving 10 years, 4.669 weeks to 10 years, 6.812 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 270/500:\x1B[0m Solving 10 years, 6.812 weeks to 10 years, 8.955 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 271/500:\x1B[0m Solving 10 years, 8.955 weeks to 10 years, 11.1 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 272/500:\x1B[0m Solving 10 years, 11.1 weeks to 10 years, 13.24 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 273/500:\x1B[0m Solving 10 years, 13.24 weeks to 10 years, 15.38 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 274/500:\x1B[0m Solving 10 years, 15.38 weeks to 10 years, 17.53 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 275/500:\x1B[0m Solving 10 years, 17.53 weeks to 10 years, 19.67 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 276/500:\x1B[0m Solving 10 years, 19.67 weeks to 10 years, 21.81 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 277/500:\x1B[0m Solving 10 years, 21.81 weeks to 10 years, 23.95 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 278/500:\x1B[0m Solving 10 years, 23.95 weeks to 10 years, 26.1 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 279/500:\x1B[0m Solving 10 years, 26.1 weeks to 10 years, 28.24 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 280/500:\x1B[0m Solving 10 years, 28.24 weeks to 10 years, 30.38 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 281/500:\x1B[0m Solving 10 years, 30.38 weeks to 10 years, 32.53 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 282/500:\x1B[0m Solving 10 years, 32.53 weeks to 10 years, 34.67 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 283/500:\x1B[0m Solving 10 years, 34.67 weeks to 10 years, 36.81 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 284/500:\x1B[0m Solving 10 years, 36.81 weeks to 10 years, 38.95 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 285/500:\x1B[0m Solving 10 years, 38.95 weeks to 10 years, 41.1 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 286/500:\x1B[0m Solving 10 years, 41.1 weeks to 10 years, 43.24 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 287/500:\x1B[0m Solving 10 years, 43.24 weeks to 10 years, 45.38 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 288/500:\x1B[0m Solving 10 years, 45.38 weeks to 10 years, 47.53 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 289/500:\x1B[0m Solving 10 years, 47.53 weeks to 10 years, 49.67 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 290/500:\x1B[0m Solving 10 years, 49.67 weeks to 10 years, 51.81 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 291/500:\x1B[0m Solving 10 years, 51.81 weeks to 11 years, 1.778 week, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 292/500:\x1B[0m Solving 11 years, 1.778 week to 11 years, 3.92 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 293/500:\x1B[0m Solving 11 years, 3.92 weeks to 11 years, 6.063 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 294/500:\x1B[0m Solving 11 years, 6.063 weeks to 11 years, 8.206 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 295/500:\x1B[0m Solving 11 years, 8.206 weeks to 11 years, 10.35 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 296/500:\x1B[0m Solving 11 years, 10.35 weeks to 11 years, 12.49 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 297/500:\x1B[0m Solving 11 years, 12.49 weeks to 11 years, 14.63 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 298/500:\x1B[0m Solving 11 years, 14.63 weeks to 11 years, 16.78 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 299/500:\x1B[0m Solving 11 years, 16.78 weeks to 11 years, 18.92 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 300/500:\x1B[0m Solving 11 years, 18.92 weeks to 11 years, 21.06 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 301/500:\x1B[0m Solving 11 years, 21.06 weeks to 11 years, 23.21 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 302/500:\x1B[0m Solving 11 years, 23.21 weeks to 11 years, 25.35 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 303/500:\x1B[0m Solving 11 years, 25.35 weeks to 11 years, 27.49 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 304/500:\x1B[0m Solving 11 years, 27.49 weeks to 11 years, 29.63 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 305/500:\x1B[0m Solving 11 years, 29.63 weeks to 11 years, 31.78 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 306/500:\x1B[0m Solving 11 years, 31.78 weeks to 11 years, 33.92 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 307/500:\x1B[0m Solving 11 years, 33.92 weeks to 11 years, 36.06 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 308/500:\x1B[0m Solving 11 years, 36.06 weeks to 11 years, 38.21 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 309/500:\x1B[0m Solving 11 years, 38.21 weeks to 11 years, 40.35 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 310/500:\x1B[0m Solving 11 years, 40.35 weeks to 11 years, 42.49 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 311/500:\x1B[0m Solving 11 years, 42.49 weeks to 11 years, 44.63 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 312/500:\x1B[0m Solving 11 years, 44.63 weeks to 11 years, 46.78 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 313/500:\x1B[0m Solving 11 years, 46.78 weeks to 11 years, 48.92 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 314/500:\x1B[0m Solving 11 years, 48.92 weeks to 11 years, 51.06 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 315/500:\x1B[0m Solving 11 years, 51.06 weeks to 12 years, 1.029 week, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 316/500:\x1B[0m Solving 12 years, 1.029 week to 12 years, 3.171 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 317/500:\x1B[0m Solving 12 years, 3.171 weeks to 12 years, 5.314 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 318/500:\x1B[0m Solving 12 years, 5.314 weeks to 12 years, 7.457 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 319/500:\x1B[0m Solving 12 years, 7.457 weeks to 12 years, 9.6 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 320/500:\x1B[0m Solving 12 years, 9.6 weeks to 12 years, 11.74 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 321/500:\x1B[0m Solving 12 years, 11.74 weeks to 12 years, 13.89 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 322/500:\x1B[0m Solving 12 years, 13.89 weeks to 12 years, 16.03 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 323/500:\x1B[0m Solving 12 years, 16.03 weeks to 12 years, 18.17 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 324/500:\x1B[0m Solving 12 years, 18.17 weeks to 12 years, 20.31 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 325/500:\x1B[0m Solving 12 years, 20.31 weeks to 12 years, 22.46 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 326/500:\x1B[0m Solving 12 years, 22.46 weeks to 12 years, 24.6 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 327/500:\x1B[0m Solving 12 years, 24.6 weeks to 12 years, 26.74 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 328/500:\x1B[0m Solving 12 years, 26.74 weeks to 12 years, 28.89 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 329/500:\x1B[0m Solving 12 years, 28.89 weeks to 12 years, 31.03 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 330/500:\x1B[0m Solving 12 years, 31.03 weeks to 12 years, 33.17 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 331/500:\x1B[0m Solving 12 years, 33.17 weeks to 12 years, 35.31 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 332/500:\x1B[0m Solving 12 years, 35.31 weeks to 12 years, 37.46 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 333/500:\x1B[0m Solving 12 years, 37.46 weeks to 12 years, 39.6 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 334/500:\x1B[0m Solving 12 years, 39.6 weeks to 12 years, 41.74 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 335/500:\x1B[0m Solving 12 years, 41.74 weeks to 12 years, 43.89 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 336/500:\x1B[0m Solving 12 years, 43.89 weeks to 12 years, 46.03 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 337/500:\x1B[0m Solving 12 years, 46.03 weeks to 12 years, 48.17 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 338/500:\x1B[0m Solving 12 years, 48.17 weeks to 12 years, 50.31 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 339/500:\x1B[0m Solving 12 years, 50.31 weeks to 13 years, 1.958 day, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 340/500:\x1B[0m Solving 13 years, 1.958 day to 13 years, 2.422 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 341/500:\x1B[0m Solving 13 years, 2.422 weeks to 13 years, 4.565 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 342/500:\x1B[0m Solving 13 years, 4.565 weeks to 13 years, 6.708 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 343/500:\x1B[0m Solving 13 years, 6.708 weeks to 13 years, 8.851 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 344/500:\x1B[0m Solving 13 years, 8.851 weeks to 13 years, 10.99 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 345/500:\x1B[0m Solving 13 years, 10.99 weeks to 13 years, 13.14 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 346/500:\x1B[0m Solving 13 years, 13.14 weeks to 13 years, 15.28 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 347/500:\x1B[0m Solving 13 years, 15.28 weeks to 13 years, 17.42 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 348/500:\x1B[0m Solving 13 years, 17.42 weeks to 13 years, 19.57 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 349/500:\x1B[0m Solving 13 years, 19.57 weeks to 13 years, 21.71 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 350/500:\x1B[0m Solving 13 years, 21.71 weeks to 13 years, 23.85 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 351/500:\x1B[0m Solving 13 years, 23.85 weeks to 13 years, 25.99 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 352/500:\x1B[0m Solving 13 years, 25.99 weeks to 13 years, 28.14 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 353/500:\x1B[0m Solving 13 years, 28.14 weeks to 13 years, 30.28 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 354/500:\x1B[0m Solving 13 years, 30.28 weeks to 13 years, 32.42 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 355/500:\x1B[0m Solving 13 years, 32.42 weeks to 13 years, 34.57 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 356/500:\x1B[0m Solving 13 years, 34.57 weeks to 13 years, 36.71 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 357/500:\x1B[0m Solving 13 years, 36.71 weeks to 13 years, 38.85 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 358/500:\x1B[0m Solving 13 years, 38.85 weeks to 13 years, 40.99 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 359/500:\x1B[0m Solving 13 years, 40.99 weeks to 13 years, 43.14 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 360/500:\x1B[0m Solving 13 years, 43.14 weeks to 13 years, 45.28 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 361/500:\x1B[0m Solving 13 years, 45.28 weeks to 13 years, 47.42 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 362/500:\x1B[0m Solving 13 years, 47.42 weeks to 13 years, 49.57 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 363/500:\x1B[0m Solving 13 years, 49.57 weeks to 13 years, 51.71 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 364/500:\x1B[0m Solving 13 years, 51.71 weeks to 14 years, 1.674 week, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 365/500:\x1B[0m Solving 14 years, 1.674 week to 14 years, 3.816 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 366/500:\x1B[0m Solving 14 years, 3.816 weeks to 14 years, 5.959 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 367/500:\x1B[0m Solving 14 years, 5.959 weeks to 14 years, 8.102 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 368/500:\x1B[0m Solving 14 years, 8.102 weeks to 14 years, 10.24 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 369/500:\x1B[0m Solving 14 years, 10.24 weeks to 14 years, 12.39 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 370/500:\x1B[0m Solving 14 years, 12.39 weeks to 14 years, 14.53 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 371/500:\x1B[0m Solving 14 years, 14.53 weeks to 14 years, 16.67 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 372/500:\x1B[0m Solving 14 years, 16.67 weeks to 14 years, 18.82 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 373/500:\x1B[0m Solving 14 years, 18.82 weeks to 14 years, 20.96 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 374/500:\x1B[0m Solving 14 years, 20.96 weeks to 14 years, 23.1 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 375/500:\x1B[0m Solving 14 years, 23.1 weeks to 14 years, 25.25 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 376/500:\x1B[0m Solving 14 years, 25.25 weeks to 14 years, 27.39 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 377/500:\x1B[0m Solving 14 years, 27.39 weeks to 14 years, 29.53 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 378/500:\x1B[0m Solving 14 years, 29.53 weeks to 14 years, 31.67 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 379/500:\x1B[0m Solving 14 years, 31.67 weeks to 14 years, 33.82 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 380/500:\x1B[0m Solving 14 years, 33.82 weeks to 14 years, 35.96 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 381/500:\x1B[0m Solving 14 years, 35.96 weeks to 14 years, 38.1 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 382/500:\x1B[0m Solving 14 years, 38.1 weeks to 14 years, 40.24 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 383/500:\x1B[0m Solving 14 years, 40.24 weeks to 14 years, 42.39 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 384/500:\x1B[0m Solving 14 years, 42.39 weeks to 14 years, 44.53 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 385/500:\x1B[0m Solving 14 years, 44.53 weeks to 14 years, 46.67 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 386/500:\x1B[0m Solving 14 years, 46.67 weeks to 14 years, 48.82 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 387/500:\x1B[0m Solving 14 years, 48.82 weeks to 14 years, 50.96 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 388/500:\x1B[0m Solving 14 years, 50.96 weeks to 15 years, 6.473 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 389/500:\x1B[0m Solving 15 years, 6.473 days to 15 years, 3.067 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 390/500:\x1B[0m Solving 15 years, 3.067 weeks to 15 years, 5.21 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 391/500:\x1B[0m Solving 15 years, 5.21 weeks to 15 years, 7.353 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 392/500:\x1B[0m Solving 15 years, 7.353 weeks to 15 years, 9.496 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 393/500:\x1B[0m Solving 15 years, 9.496 weeks to 15 years, 11.64 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 394/500:\x1B[0m Solving 15 years, 11.64 weeks to 15 years, 13.78 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 395/500:\x1B[0m Solving 15 years, 13.78 weeks to 15 years, 15.92 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 396/500:\x1B[0m Solving 15 years, 15.92 weeks to 15 years, 18.07 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 397/500:\x1B[0m Solving 15 years, 18.07 weeks to 15 years, 20.21 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 398/500:\x1B[0m Solving 15 years, 20.21 weeks to 15 years, 22.35 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 399/500:\x1B[0m Solving 15 years, 22.35 weeks to 15 years, 24.5 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 400/500:\x1B[0m Solving 15 years, 24.5 weeks to 15 years, 26.64 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 401/500:\x1B[0m Solving 15 years, 26.64 weeks to 15 years, 28.78 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 402/500:\x1B[0m Solving 15 years, 28.78 weeks to 15 years, 30.92 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 403/500:\x1B[0m Solving 15 years, 30.92 weeks to 15 years, 33.07 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 404/500:\x1B[0m Solving 15 years, 33.07 weeks to 15 years, 35.21 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 405/500:\x1B[0m Solving 15 years, 35.21 weeks to 15 years, 37.35 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 406/500:\x1B[0m Solving 15 years, 37.35 weeks to 15 years, 39.5 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 407/500:\x1B[0m Solving 15 years, 39.5 weeks to 15 years, 41.64 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 408/500:\x1B[0m Solving 15 years, 41.64 weeks to 15 years, 43.78 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 409/500:\x1B[0m Solving 15 years, 43.78 weeks to 15 years, 45.92 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 410/500:\x1B[0m Solving 15 years, 45.92 weeks to 15 years, 48.07 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 411/500:\x1B[0m Solving 15 years, 48.07 weeks to 15 years, 50.21 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 412/500:\x1B[0m Solving 15 years, 50.21 weeks to 16 years, 1.23 day, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 413/500:\x1B[0m Solving 16 years, 1.23 day to 16 years, 2.319 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 414/500:\x1B[0m Solving 16 years, 2.319 weeks to 16 years, 4.461 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 415/500:\x1B[0m Solving 16 years, 4.461 weeks to 16 years, 6.604 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 416/500:\x1B[0m Solving 16 years, 6.604 weeks to 16 years, 8.747 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 417/500:\x1B[0m Solving 16 years, 8.747 weeks to 16 years, 10.89 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 418/500:\x1B[0m Solving 16 years, 10.89 weeks to 16 years, 13.03 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 419/500:\x1B[0m Solving 16 years, 13.03 weeks to 16 years, 15.18 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 420/500:\x1B[0m Solving 16 years, 15.18 weeks to 16 years, 17.32 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 421/500:\x1B[0m Solving 16 years, 17.32 weeks to 16 years, 19.46 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 422/500:\x1B[0m Solving 16 years, 19.46 weeks to 16 years, 21.6 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 423/500:\x1B[0m Solving 16 years, 21.6 weeks to 16 years, 23.75 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 424/500:\x1B[0m Solving 16 years, 23.75 weeks to 16 years, 25.89 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 425/500:\x1B[0m Solving 16 years, 25.89 weeks to 16 years, 28.03 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 426/500:\x1B[0m Solving 16 years, 28.03 weeks to 16 years, 30.18 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 427/500:\x1B[0m Solving 16 years, 30.18 weeks to 16 years, 32.32 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 428/500:\x1B[0m Solving 16 years, 32.32 weeks to 16 years, 34.46 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 429/500:\x1B[0m Solving 16 years, 34.46 weeks to 16 years, 36.6 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 430/500:\x1B[0m Solving 16 years, 36.6 weeks to 16 years, 38.75 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 431/500:\x1B[0m Solving 16 years, 38.75 weeks to 16 years, 40.89 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 432/500:\x1B[0m Solving 16 years, 40.89 weeks to 16 years, 43.03 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 433/500:\x1B[0m Solving 16 years, 43.03 weeks to 16 years, 45.18 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 434/500:\x1B[0m Solving 16 years, 45.18 weeks to 16 years, 47.32 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 435/500:\x1B[0m Solving 16 years, 47.32 weeks to 16 years, 49.46 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 436/500:\x1B[0m Solving 16 years, 49.46 weeks to 16 years, 51.6 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 437/500:\x1B[0m Solving 16 years, 51.6 weeks to 17 years, 1.57 week, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 438/500:\x1B[0m Solving 17 years, 1.57 week to 17 years, 3.712 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 439/500:\x1B[0m Solving 17 years, 3.712 weeks to 17 years, 5.855 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 440/500:\x1B[0m Solving 17 years, 5.855 weeks to 17 years, 7.998 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 441/500:\x1B[0m Solving 17 years, 7.998 weeks to 17 years, 10.14 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 442/500:\x1B[0m Solving 17 years, 10.14 weeks to 17 years, 12.28 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 443/500:\x1B[0m Solving 17 years, 12.28 weeks to 17 years, 14.43 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 444/500:\x1B[0m Solving 17 years, 14.43 weeks to 17 years, 16.57 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 445/500:\x1B[0m Solving 17 years, 16.57 weeks to 17 years, 18.71 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 446/500:\x1B[0m Solving 17 years, 18.71 weeks to 17 years, 20.86 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 447/500:\x1B[0m Solving 17 years, 20.86 weeks to 17 years, 23 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 448/500:\x1B[0m Solving 17 years, 23 weeks to 17 years, 25.14 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 449/500:\x1B[0m Solving 17 years, 25.14 weeks to 17 years, 27.28 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 450/500:\x1B[0m Solving 17 years, 27.28 weeks to 17 years, 29.43 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 451/500:\x1B[0m Solving 17 years, 29.43 weeks to 17 years, 31.57 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 452/500:\x1B[0m Solving 17 years, 31.57 weeks to 17 years, 33.71 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 453/500:\x1B[0m Solving 17 years, 33.71 weeks to 17 years, 35.86 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 454/500:\x1B[0m Solving 17 years, 35.86 weeks to 17 years, 38 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 455/500:\x1B[0m Solving 17 years, 38 weeks to 17 years, 40.14 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 456/500:\x1B[0m Solving 17 years, 40.14 weeks to 17 years, 42.28 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 457/500:\x1B[0m Solving 17 years, 42.28 weeks to 17 years, 44.43 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 458/500:\x1B[0m Solving 17 years, 44.43 weeks to 17 years, 46.57 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 459/500:\x1B[0m Solving 17 years, 46.57 weeks to 17 years, 48.71 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 460/500:\x1B[0m Solving 17 years, 48.71 weeks to 17 years, 50.86 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 461/500:\x1B[0m Solving 17 years, 50.86 weeks to 18 years, 5.745 days, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 462/500:\x1B[0m Solving 18 years, 5.745 days to 18 years, 2.964 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 463/500:\x1B[0m Solving 18 years, 2.964 weeks to 18 years, 5.106 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 464/500:\x1B[0m Solving 18 years, 5.106 weeks to 18 years, 7.249 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 465/500:\x1B[0m Solving 18 years, 7.249 weeks to 18 years, 9.392 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 466/500:\x1B[0m Solving 18 years, 9.392 weeks to 18 years, 11.54 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 467/500:\x1B[0m Solving 18 years, 11.54 weeks to 18 years, 13.68 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 468/500:\x1B[0m Solving 18 years, 13.68 weeks to 18 years, 15.82 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 469/500:\x1B[0m Solving 18 years, 15.82 weeks to 18 years, 17.96 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 470/500:\x1B[0m Solving 18 years, 17.96 weeks to 18 years, 20.11 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 471/500:\x1B[0m Solving 18 years, 20.11 weeks to 18 years, 22.25 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 472/500:\x1B[0m Solving 18 years, 22.25 weeks to 18 years, 24.39 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 473/500:\x1B[0m Solving 18 years, 24.39 weeks to 18 years, 26.54 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 474/500:\x1B[0m Solving 18 years, 26.54 weeks to 18 years, 28.68 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 475/500:\x1B[0m Solving 18 years, 28.68 weeks to 18 years, 30.82 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 476/500:\x1B[0m Solving 18 years, 30.82 weeks to 18 years, 32.96 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 477/500:\x1B[0m Solving 18 years, 32.96 weeks to 18 years, 35.11 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 478/500:\x1B[0m Solving 18 years, 35.11 weeks to 18 years, 37.25 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 479/500:\x1B[0m Solving 18 years, 37.25 weeks to 18 years, 39.39 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 480/500:\x1B[0m Solving 18 years, 39.39 weeks to 18 years, 41.53 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 481/500:\x1B[0m Solving 18 years, 41.53 weeks to 18 years, 43.68 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 482/500:\x1B[0m Solving 18 years, 43.68 weeks to 18 years, 45.82 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 483/500:\x1B[0m Solving 18 years, 45.82 weeks to 18 years, 47.96 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 484/500:\x1B[0m Solving 18 years, 47.96 weeks to 18 years, 50.11 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 485/500:\x1B[0m Solving 18 years, 50.11 weeks to 19 years, 12.06 hours, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 486/500:\x1B[0m Solving 19 years, 12.06 hours to 19 years, 2.215 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 487/500:\x1B[0m Solving 19 years, 2.215 weeks to 19 years, 4.357 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 488/500:\x1B[0m Solving 19 years, 4.357 weeks to 19 years, 6.5 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 489/500:\x1B[0m Solving 19 years, 6.5 weeks to 19 years, 8.643 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 490/500:\x1B[0m Solving 19 years, 8.643 weeks to 19 years, 10.79 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 491/500:\x1B[0m Solving 19 years, 10.79 weeks to 19 years, 12.93 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 492/500:\x1B[0m Solving 19 years, 12.93 weeks to 19 years, 15.07 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 493/500:\x1B[0m Solving 19 years, 15.07 weeks to 19 years, 17.21 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 494/500:\x1B[0m Solving 19 years, 17.21 weeks to 19 years, 19.36 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 495/500:\x1B[0m Solving 19 years, 19.36 weeks to 19 years, 21.5 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 496/500:\x1B[0m Solving 19 years, 21.5 weeks to 19 years, 23.64 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 497/500:\x1B[0m Solving 19 years, 23.64 weeks to 19 years, 25.79 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 498/500:\x1B[0m Solving 19 years, 25.79 weeks to 19 years, 27.93 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 499/500:\x1B[0m Solving 19 years, 27.93 weeks to 19 years, 30.07 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[34;1mStep 500/500:\x1B[0m Solving 19 years, 30.07 weeks to 19 years, 32.21 weeks, Δt = 2 weeks, 1 day</span></span>
<span class="line"><span>\x1B[92;1;4mSimulation complete:\x1B[0m Completed 500 report steps in 35 seconds, 130 milliseconds, 661 microseconds and 1526 iterations.</span></span>
<span class="line"><span>╭────────────────┬───────────┬───────────────┬──────────╮</span></span>
<span class="line"><span>│ Iteration type │  Avg/step │  Avg/ministep │    Total │</span></span>
<span class="line"><span>│                │ 500 steps │ 503 ministeps │ (wasted) │</span></span>
<span class="line"><span>├────────────────┼───────────┼───────────────┼──────────┤</span></span>
<span class="line"><span>│ Newton         │     3.052 │        3.0338 │ 1526 (0) │</span></span>
<span class="line"><span>│ Linearization  │     4.058 │        4.0338 │ 2029 (0) │</span></span>
<span class="line"><span>│ Linear solver  │     3.052 │        3.0338 │ 1526 (0) │</span></span>
<span class="line"><span>│ Precond apply  │     6.104 │       6.06759 │ 3052 (0) │</span></span>
<span class="line"><span>╰────────────────┴───────────┴───────────────┴──────────╯</span></span>
<span class="line"><span>╭───────────────┬─────────┬────────────┬─────────╮</span></span>
<span class="line"><span>│ Timing type   │    Each │   Relative │   Total │</span></span>
<span class="line"><span>│               │      ms │ Percentage │       s │</span></span>
<span class="line"><span>├───────────────┼─────────┼────────────┼─────────┤</span></span>
<span class="line"><span>│ Properties    │ 17.9467 │    77.96 % │ 27.3866 │</span></span>
<span class="line"><span>│ Equations     │  1.0096 │     5.83 % │  2.0485 │</span></span>
<span class="line"><span>│ Assembly      │  0.5545 │     3.20 % │  1.1251 │</span></span>
<span class="line"><span>│ Linear solve  │  0.2060 │     0.89 % │  0.3144 │</span></span>
<span class="line"><span>│ Linear setup  │  0.7144 │     3.10 % │  1.0902 │</span></span>
<span class="line"><span>│ Precond apply │  0.1069 │     0.93 % │  0.3262 │</span></span>
<span class="line"><span>│ Update        │  0.1877 │     0.82 % │  0.2865 │</span></span>
<span class="line"><span>│ Convergence   │  0.6080 │     3.51 % │  1.2336 │</span></span>
<span class="line"><span>│ Input/Output  │  0.3337 │     0.48 % │  0.1679 │</span></span>
<span class="line"><span>│ Other         │  0.7547 │     3.28 % │  1.1517 │</span></span>
<span class="line"><span>├───────────────┼─────────┼────────────┼─────────┤</span></span>
<span class="line"><span>│ Total         │ 23.0214 │   100.00 % │ 35.1307 │</span></span>
<span class="line"><span>╰───────────────┴─────────┴────────────┴─────────╯</span></span></code></pre></div><h2 id="Plot-solutions-and-compare" tabindex="-1">Plot solutions and compare <a class="header-anchor" href="#Plot-solutions-and-compare" aria-label="Permalink to &quot;Plot solutions and compare {#Plot-solutions-and-compare}&quot;">​</a></h2><p>The 1D displacement can be plotted as a line plot. We pick a step midway through the simulation and plot compositions, saturations and pressure.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">cmap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :tableau_hue_circle</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ref_path </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> joinpath</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(dpth, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;reference.jld2&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ref </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Jutul</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">JLD2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">load</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ref_path)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">step_to_plot </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 250</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> with_theme</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">theme_latexfonts</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">do</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reservoir_domain</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case)[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:cell_centroids</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    mz </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 3</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ix </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> step_to_plot</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    mt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :circle</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">800</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">400</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], xlabel </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Cell center / m&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cnames </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;DECANE&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;CO2&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;METHANE&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cnames </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;C₁₀&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;CO₂&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;C₁&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cnames </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\t</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ext{C}_{10}&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\t</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ext{CO}_2&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\t</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ext{C}_1&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    lineh </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> []</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    lnames </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> []</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    crange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> range</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(crange</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 4</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            cname </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\t</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ext{S}_g&quot;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            gprs </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> missing</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            ecl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ref[</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;e300&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][ix][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Saturations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            ju </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> states[ix][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Saturations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        else</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">            @assert</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 4</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            ecl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ref[</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;e300&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][ix][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:OverallMoleFractions</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][i, :]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            gprs </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ref[</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;adgprs&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][ix][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:OverallMoleFractions</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][i, :]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            ju </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> states[ix][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:OverallMoleFractions</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][i, :]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            cname </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cnames[i]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        h </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, x, ju, colormap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cmap, color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> crange, label </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cname)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        push!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(lnames, cname)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        push!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(lineh, h)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        scatter!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, x, ecl, markersize </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> mz, colormap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cmap, color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> crange)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        if</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> !</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ismissing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(gprs)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">            lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, x, gprs, colormap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cmap, color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i, colorrange </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> crange, linestyle </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :dash</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    l_ju </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> LineElement</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :black</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, linestyle </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> nothing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    l_ecl </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> MarkerElement</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :black</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, markersize </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> mz, marker </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> mt)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    l_gprs </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> LineElement</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :black</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, linestyle </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :dash</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    Legend</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">],</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        [[l_ju, l_ecl, l_gprs], lineh],</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        [[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\t</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ext{JutulDarcy}&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\t</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ext{E300}&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\t</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ext{AD-GPRS}&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], lnames],</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        [</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Simulator&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Result&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">],</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        tellwidth </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        orientation </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :horizontal</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    )</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    fig</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><p><img src="`+l+`" alt=""></p><h2 id="Calculate-sensitivities" tabindex="-1">Calculate sensitivities <a class="header-anchor" href="#Calculate-sensitivities" aria-label="Permalink to &quot;Calculate sensitivities {#Calculate-sensitivities}&quot;">​</a></h2><p>We demonstrate how the parameter sensitivities of an objective function can be calculated for a compositional model.</p><p>The objective function is taken to be the average gas saturation at a specific report step that was plotted in the previous section.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">import</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Statistics</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> mean</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">import</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> JutulDarcy</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> reservoir_sensitivities</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> objective_function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(model, state, Δt, step_i, forces)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> step_i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">!=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> step_to_plot</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        return</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.0</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    sg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> @view</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> state</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Reservoir</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Saturations[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> mean</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sg)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">data_domain_with_gradients </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reservoir_sensitivities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case, result, objective_function, include_parameters </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> with_theme</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">theme_latexfonts</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">do</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reservoir_domain</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(case)[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:cell_centroids</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">][</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, :]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    mz </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 3</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ix </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> step_to_plot</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cmap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Dark2_5</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cmap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Accent_4</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cmap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :Spectral_4</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cmap </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :tableau_hue_circle</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    mt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :circle</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">800</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">400</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    normalize </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">./</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">maximum</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> minimum</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    logscale </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sign</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(x)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">log10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">abs</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(x))</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ∂T </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data_domain_with_gradients[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:temperature</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ∂ϕ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> data_domain_with_gradients[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:porosity</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ax1 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], yticklabelcolor </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :blue</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, xlabel </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Cell center / m&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ax2 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], yticklabelcolor </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :red</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, yaxisposition </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :right</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    hidespines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax2)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    hidexdecorations!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax2)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    l1 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax1, x, ∂T, label </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Temperature&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :blue</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    l2 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax2, x, ∂ϕ, label </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Porosity&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :red</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    Legend</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], [l1, l2], [</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Temperature&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Porosity&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Parameter&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, tellwidth </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, orientation </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :horizontal</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    fig</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+k+'" alt=""></p><h2 id="Example-on-GitHub" tabindex="-1">Example on GitHub <a class="header-anchor" href="#Example-on-GitHub" aria-label="Permalink to &quot;Example on GitHub {#Example-on-GitHub}&quot;">​</a></h2><p>If you would like to run this example yourself, it can be downloaded from the JutulDarcy.jl GitHub repository <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/main/examples/validation_compositional.jl" target="_blank" rel="noreferrer">as a script</a>, or as a <a href="https://github.com/sintefmath/JutulDarcy.jl/blob/gh-pages/dev/final_site/notebooks/validation_compositional.ipynb" target="_blank" rel="noreferrer">Jupyter Notebook</a></p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',15))])}const v=t(y,[["render",u]]);export{c as __pageData,v as default};
