import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/JutulDarcy.jl/dev/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Manual', collapsed: false, items: [
{ text: 'Introduction', collapsed: false, items: [
{ text: 'JutulDarcy.jl', link: '/index' },
{ text: 'Getting started', link: '/man/intro' },
{ text: 'Your first JutulDarcy.jl simulation', link: '/man/first_ex' },
{ text: 'FAQ', link: '/extras/faq' }]
 },
{ text: 'Fundamentals', collapsed: false, items: [
{ text: 'High-level API', link: '/man/highlevel' },
{ text: 'Input formats', link: '/man/basics/input_files' },
{ text: 'Supported physical systems', link: '/man/basics/systems' },
{ text: 'Solving the equations', link: '/man/basics/solution' }]
 },
{ text: 'Detailed API', collapsed: false, items: [
{ text: 'Driving forces', link: '/man/basics/forces' },
{ text: 'Wells and controls', link: '/man/basics/wells' },
{ text: 'Primary variables', link: '/man/basics/primary' },
{ text: 'Secondary variables (properties)', link: '/man/basics/secondary' },
{ text: 'Parameters', link: '/man/basics/parameters' },
{ text: 'Plotting and visualization', link: '/man/basics/plotting' },
{ text: 'Utilities', link: '/man/basics/utilities' }]
 },
{ text: 'Parallelism and compilation', collapsed: false, items: [
{ text: 'Multi-threading and MPI support', link: '/man/advanced/mpi' },
{ text: 'GPU support', link: '/man/advanced/gpu' },
{ text: 'Standalone reservoir simulator', link: '/man/advanced/compiled' }]
 },
{ text: 'References', collapsed: false, items: [
{ text: 'Package docstring', link: '/man/basics/package' },
{ text: 'Jutul functions', link: '/ref/jutul' },
{ text: 'Bibliography', link: '/extras/refs' }]
 }]
 },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Introduction', collapsed: false, items: [
{ text: 'Simulating Eclipse/DATA input files', link: '/examples/introduction/data_input_file' },
{ text: 'Intro to sensitivities in JutulDarcy', link: '/examples/introduction/intro_sensitivities' },
{ text: 'Buckley-Leverett two-phase problem', link: '/examples/introduction/two_phase_buckley_leverett' },
{ text: 'Gravity segregation example', link: '/examples/introduction/two_phase_gravity_segregation' },
{ text: 'Gravity circulation with CPR preconditioner', link: '/examples/introduction/two_phase_unstable_gravity' },
{ text: 'Introduction to wells', link: '/examples/introduction/wells_intro' }]
 },
{ text: 'Workflow', collapsed: false, items: [
{ text: 'Adding new wells to an existing model', link: '/examples/workflow/adding_new_wells' },
{ text: 'CO2 injection in saline aquifer with storage inventory', link: '/examples/workflow/co2_sloped' },
{ text: 'Hydrostatic equilibriation of models', link: '/examples/workflow/equilibrium_state' },
{ text: 'Quarter-five-spot example', link: '/examples/workflow/five_spot_ensemble' },
{ text: 'Hybrid simulation with neural network for relative permeability', link: '/examples/workflow/hybrid_simulation_relperm' },
{ text: 'Model coarsening', link: '/examples/workflow/model_coarsening' },
{ text: 'Gradient-based optimization of net present value (NPV)', link: '/examples/workflow/rate_optimization' },
{ text: 'Adding tracers to a flow simulation', link: '/examples/workflow/tracers_two_wells' }]
 },
{ text: 'Data assimilation', collapsed: false, items: [
{ text: 'Advanced history matching: Regions, blending and parametric models', link: '/examples/data_assimilation/advanced_history_match' },
{ text: 'History matching a coarse model - CGNet', link: '/examples/data_assimilation/cgnet_egg' },
{ text: 'Gradient-based matching of parameters against observations', link: '/examples/data_assimilation/optimize_simple_bl' }]
 },
{ text: 'Geothermal', collapsed: false, items: [
{ text: 'Borehole Thermal Energy Storage (BTES)', link: '/examples/geothermal/btes' },
{ text: 'Geothermal doublet', link: '/examples/geothermal/geothermal_doublet' },
{ text: 'High-temperature Aquifer Thermal Energy Storage (HT-ATES)', link: '/examples/geothermal/htates_intro' }]
 },
{ text: 'Compositional', collapsed: false, items: [
{ text: 'Intro to compositional flow', link: '/examples/compositional/compositional_2d_vertical' },
{ text: 'A more complex compositional model', link: '/examples/compositional/compositional_5components' }]
 },
{ text: 'Discretization', collapsed: false, items: [
{ text: 'Consistent discretizations: Average MPFA and nonlinear TPFA', link: '/examples/discretization/consistent_avgmpfa' },
{ text: 'Consistent and high-resolution: WENO, NTPFA and AvgMPFA', link: '/examples/discretization/mpfa_weno_discretizations' }]
 },
{ text: 'Properties', collapsed: false, items: [
{ text: 'CO2-brine correlations with salinity', link: '/examples/properties/co2_props' },
{ text: 'Relative Permeabilities in JutulDarcy', link: '/examples/properties/relperms' }]
 }]
 },
{ text: 'Validation', collapsed: false, items: [
{ text: 'Overview', link: '/man/validation' },
{ text: 'Models', collapsed: false, items: [
{ text: 'Validation of equation-of-state compositional flow', link: '/examples/validation/validation_compositional' },
{ text: 'The Egg model: Two-phase oil-water model', link: '/examples/validation/validation_egg' },
{ text: 'Comparison between JutulDarcy.jl and MRST', link: '/examples/validation/validation_mrst' },
{ text: 'Norne: Real field black-oil model', link: '/examples/validation/validation_norne_nohyst' },
{ text: 'The OLYMPUS benchmark model: Two-phase corner-point reservoir', link: '/examples/validation/validation_olympus_1' },
{ text: 'Polymer injection in a 2D black-oil reservoir model', link: '/examples/validation/validation_polymer' },
{ text: 'SPE1: Small black-oil gas injection', link: '/examples/validation/validation_spe1' },
{ text: 'SPE9: Black-oil depletion with dissolved gas', link: '/examples/validation/validation_spe9' },
{ text: 'Aquifer thermal energy storage (ATES) validation', link: '/examples/validation/validation_thermal' }]
 }]
 }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/JutulDarcy.jl/dev/',// TODO: replace this in makedocs!
  title: 'JutulDarcy.jl',
  description: 'Documentation for JutulDarcy.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: '/favicon.ico' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  ignoreDeadLinks: true,
  vite: {
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    logo: { src: '/logo.png', width: 24, height: 24},
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Manual', collapsed: false, items: [
{ text: 'Introduction', collapsed: false, items: [
{ text: 'JutulDarcy.jl', link: '/index' },
{ text: 'Getting started', link: '/man/intro' },
{ text: 'Your first JutulDarcy.jl simulation', link: '/man/first_ex' },
{ text: 'FAQ', link: '/extras/faq' }]
 },
{ text: 'Fundamentals', collapsed: false, items: [
{ text: 'High-level API', link: '/man/highlevel' },
{ text: 'Input formats', link: '/man/basics/input_files' },
{ text: 'Supported physical systems', link: '/man/basics/systems' },
{ text: 'Solving the equations', link: '/man/basics/solution' }]
 },
{ text: 'Detailed API', collapsed: false, items: [
{ text: 'Driving forces', link: '/man/basics/forces' },
{ text: 'Wells and controls', link: '/man/basics/wells' },
{ text: 'Primary variables', link: '/man/basics/primary' },
{ text: 'Secondary variables (properties)', link: '/man/basics/secondary' },
{ text: 'Parameters', link: '/man/basics/parameters' },
{ text: 'Plotting and visualization', link: '/man/basics/plotting' },
{ text: 'Utilities', link: '/man/basics/utilities' }]
 },
{ text: 'Parallelism and compilation', collapsed: false, items: [
{ text: 'Multi-threading and MPI support', link: '/man/advanced/mpi' },
{ text: 'GPU support', link: '/man/advanced/gpu' },
{ text: 'Standalone reservoir simulator', link: '/man/advanced/compiled' }]
 },
{ text: 'References', collapsed: false, items: [
{ text: 'Package docstring', link: '/man/basics/package' },
{ text: 'Jutul functions', link: '/ref/jutul' },
{ text: 'Bibliography', link: '/extras/refs' }]
 }]
 },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Introduction', collapsed: false, items: [
{ text: 'Simulating Eclipse/DATA input files', link: '/examples/introduction/data_input_file' },
{ text: 'Intro to sensitivities in JutulDarcy', link: '/examples/introduction/intro_sensitivities' },
{ text: 'Buckley-Leverett two-phase problem', link: '/examples/introduction/two_phase_buckley_leverett' },
{ text: 'Gravity segregation example', link: '/examples/introduction/two_phase_gravity_segregation' },
{ text: 'Gravity circulation with CPR preconditioner', link: '/examples/introduction/two_phase_unstable_gravity' },
{ text: 'Introduction to wells', link: '/examples/introduction/wells_intro' }]
 },
{ text: 'Workflow', collapsed: false, items: [
{ text: 'Adding new wells to an existing model', link: '/examples/workflow/adding_new_wells' },
{ text: 'CO2 injection in saline aquifer with storage inventory', link: '/examples/workflow/co2_sloped' },
{ text: 'Hydrostatic equilibriation of models', link: '/examples/workflow/equilibrium_state' },
{ text: 'Quarter-five-spot example', link: '/examples/workflow/five_spot_ensemble' },
{ text: 'Hybrid simulation with neural network for relative permeability', link: '/examples/workflow/hybrid_simulation_relperm' },
{ text: 'Model coarsening', link: '/examples/workflow/model_coarsening' },
{ text: 'Gradient-based optimization of net present value (NPV)', link: '/examples/workflow/rate_optimization' },
{ text: 'Adding tracers to a flow simulation', link: '/examples/workflow/tracers_two_wells' }]
 },
{ text: 'Data assimilation', collapsed: false, items: [
{ text: 'Advanced history matching: Regions, blending and parametric models', link: '/examples/data_assimilation/advanced_history_match' },
{ text: 'History matching a coarse model - CGNet', link: '/examples/data_assimilation/cgnet_egg' },
{ text: 'Gradient-based matching of parameters against observations', link: '/examples/data_assimilation/optimize_simple_bl' }]
 },
{ text: 'Geothermal', collapsed: false, items: [
{ text: 'Borehole Thermal Energy Storage (BTES)', link: '/examples/geothermal/btes' },
{ text: 'Geothermal doublet', link: '/examples/geothermal/geothermal_doublet' },
{ text: 'High-temperature Aquifer Thermal Energy Storage (HT-ATES)', link: '/examples/geothermal/htates_intro' }]
 },
{ text: 'Compositional', collapsed: false, items: [
{ text: 'Intro to compositional flow', link: '/examples/compositional/compositional_2d_vertical' },
{ text: 'A more complex compositional model', link: '/examples/compositional/compositional_5components' }]
 },
{ text: 'Discretization', collapsed: false, items: [
{ text: 'Consistent discretizations: Average MPFA and nonlinear TPFA', link: '/examples/discretization/consistent_avgmpfa' },
{ text: 'Consistent and high-resolution: WENO, NTPFA and AvgMPFA', link: '/examples/discretization/mpfa_weno_discretizations' }]
 },
{ text: 'Properties', collapsed: false, items: [
{ text: 'CO2-brine correlations with salinity', link: '/examples/properties/co2_props' },
{ text: 'Relative Permeabilities in JutulDarcy', link: '/examples/properties/relperms' }]
 }]
 },
{ text: 'Validation', collapsed: false, items: [
{ text: 'Overview', link: '/man/validation' },
{ text: 'Models', collapsed: false, items: [
{ text: 'Validation of equation-of-state compositional flow', link: '/examples/validation/validation_compositional' },
{ text: 'The Egg model: Two-phase oil-water model', link: '/examples/validation/validation_egg' },
{ text: 'Comparison between JutulDarcy.jl and MRST', link: '/examples/validation/validation_mrst' },
{ text: 'Norne: Real field black-oil model', link: '/examples/validation/validation_norne_nohyst' },
{ text: 'The OLYMPUS benchmark model: Two-phase corner-point reservoir', link: '/examples/validation/validation_olympus_1' },
{ text: 'Polymer injection in a 2D black-oil reservoir model', link: '/examples/validation/validation_polymer' },
{ text: 'SPE1: Small black-oil gas injection', link: '/examples/validation/validation_spe1' },
{ text: 'SPE9: Black-oil depletion with dissolved gas', link: '/examples/validation/validation_spe9' },
{ text: 'Aquifer thermal energy storage (ATES) validation', link: '/examples/validation/validation_thermal' }]
 }]
 }
]
,
    editLink: { pattern: "https://https://github.com/sintefmath/JutulDarcy.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/sintefmath/JutulDarcy.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()} SINTEF Digital.`
    }
  }
})
