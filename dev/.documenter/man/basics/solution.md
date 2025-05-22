
# Solving the equations {#Solving-the-equations}

By default, Jutul solves a system as a fully-coupled implicit system of equations discretized with a two-point flux approximation with single-point upwind.

## Newton&#39;s method {#Newton's-method}

The standard way of solving a system of non-linear equations is by Newton&#39;s method (also known as Newton-Raphson&#39;s method). A quick recap: For a vector valued residual $\mathbf{r}(x)$ of the primary variable vector $\mathbf{x}$ we can defined a Newton update:

$$\mathbf{x}^{k+1} = \mathbf{r}^k - J^{-1} \mathbf{r}(\mathbf{x}^k), \quad J_{ij} = \frac{\partial \mathbf{r}_i}{\partial \mathbf{x}_j}.$$

JutulDarcy solves systems that generally have both non-smooth behavior and physical constraints on the values for $\textbf{x}$. For that reason, we modify Newton&#39;s method slightly:

$$\mathbf{x}^{k+1} = \mathbf{r}^k + \omega(\Delta \mathbf{x})$$

Here, $\omega$ is a function that limits the variables so that they do not change too much (e.g. Appleyard chopping, limiting of pressure, saturation and composition updates) and that they are within the prescribed limits. There are also options for automated global dampening in the presence of convergence issues. The update is then defined from inverting the Jacobian:

$$\Delta \mathbf{x} = -J^{-1} \mathbf{r}(\mathbf{x}^k), \quad J_{ij} = \frac{\partial \mathbf{r}_i}{\partial \mathbf{x}_j}.$$

Starting with $\mathbf{x}^0$as some initial guess taken from the previous time-step, we can solve the system by iterating upon this loop.

## Linear solvers and linear systems {#Linear-solvers-and-linear-systems}

For most practical applications it is not feasible or efficient to invert the Jacobian. JutulDarcy uses preconditioned iterative solvers by default, but it is possible to use direct solvers as well when working with smaller models. The high level interface for setting up a reservoir model [`setup_reservoir_model`](/man/highlevel#JutulDarcy.setup_reservoir_model) has an optional `block_backend=true` keyword argument that determines the matrix format, and consequently the linear solver type to be used.

### Direct solvers {#Direct-solvers}

If `block_backend` is set to `false`, Jutul will assemble into the standard Julia CSC sparse matrix with `Float64` elements and Julia&#39;s default direct solver will be used. It is also possible to use other Julia solvers on this system, but the default preconditioners assume that block backend is enabled.

### Iterative solver {#Iterative-solver}

If `block_backend` is set to `true`, Jutul will by default use a constrained-pressure residual (CPR) preconditioner for BiCGStab. Jutul relies on [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) for iterative solvers. The main function that selects the linear solver is [`reservoir_linsolve`](/man/basics/solution#JutulDarcy.reservoir_linsolve) that allows for the selection of different preconditioners and linear solvers. This is often an instance of [`Jutul.GenericKrylov`](/ref/jutul#Jutul.GenericKrylov) with the approprioate preconditioner.
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.reservoir_linsolve' href='#JutulDarcy.reservoir_linsolve'><span class="jlbinding">JutulDarcy.reservoir_linsolve</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
reservoir_linsolve(model, precond = :cpr; <keyword arguments>)
```


Set up iterative linear solver for a reservoir model from [`setup_reservoir_model`](/man/highlevel#JutulDarcy.setup_reservoir_model).

**Arguments**
- `model`: Reservoir model that will linearize the equations for the linear solver
  
- `precond=:cpr`: Preconditioner type to use: Either :cpr (Constrained-Pressure-Residual) or :ilu0 (block-incomplete-LU) (no effect if `solver = :direct`).
  
- `v=0`: verbosity (can lead to a large amount of output)
  
- `solver=:bicgstab`: the symbol of a Krylov.jl solver (typically :gmres or :bicgstab)
  
- `update_interval=:once`: how often the CPR AMG hierarchy is reconstructed (:once, :iteration, :ministep, :step)
  
- `update_interval_partial=:iteration`: how often the pressure system is updated in CPR
  
- `max_coarse`: max size of coarse level if using AMG
  
- `cpr_type=nothing`: type of CPR (`:true_impes`, `:quasi_impes` or `nothing` for automatic)
  
- `partial_update=true`: perform partial update of CPR preconditioner outside of AMG update (see above)
  
- `rtol=1e-3`: relative tolerance for the linear solver
  
- `max_iterations=100`: limit for linear solver iterations
  

Additional keywords are passed onto the linear solver constructor.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/linsolve.jl#L1-L20" target="_blank" rel="noreferrer">source</a></Badge>

</details>


#### Single model (only porous medium) {#Single-model-only-porous-medium}

If the model is a single model (e.g. only a reservoir) the matrix format is a block-CSC matrix that combines Julia&#39;s builtin sparse matrix format with statically sized elements from the [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) package. If we consider the two-phase immiscible system from [Multi-phase, immiscible flow](/man/basics/systems#Multi-phase,-immiscible-flow) we have a pair of equations $R_n, R_w$ together with the corresponding primary variables pressure and first saturation $p, S_n$ defined for all $N_c$ cells. Let us simplify the notation a bit so that the subscripts of the primary variables are $p, s$ and define a $N_c \times N_c$ block Jacobian linear system where the entires are given by:

$$J_{ij} = \begin{bmatrix}
   \left(\frac{\partial r_{n}}{\partial p}\right)_{ij} & \left(\frac{\partial r_{n}}{\partial s}\right)_{ij} \\
   \left(\frac{\partial r_{w}}{\partial p}\right)_{ij} & \left(\frac{\partial r_{w}}{\partial s}\right)_{ij} \end{bmatrix} = \begin{bmatrix}
   J_{np} & J_{ns} \\
   J_{wp} & J_{ws}
\end{bmatrix}_{ij}$$

This block system has several advantages:
- We immediately get access to more powerful version of standard Julia  preconditioners provided that all operations used are applicable for matrices  and are applied in the right commutative order. For example, JutulDarcy uses  the [ILUZero.jl](https://github.com/mcovalt/ILUZero.jl) package when a CSC  linear system is preconditioned with incomplete LU factorization with zero  fill-in.
  
- Sparse matrix vector products are much more efficient as less indicies need to  be looked up for each element wise multiplication.
  
- Performing local reductions over variables is much easier when they are  located in a local matrix.
  

##### Constrained Pressure Residual {#Constrained-Pressure-Residual}

The CPR preconditioner [[3](/extras/refs#wallis-cpr), [4](/extras/refs#cao-cpr)] [`CPRPreconditioner`](/man/basics/solution#JutulDarcy.CPRPreconditioner) is a multi-stage physics-informed preconditioner that seeks to decouple the global pressure part of the system from the local  transport part. In the limits of incompressible flow without gravity it can be thought of as an elliptic / hyperbolic splitting. We also implement a special variant for the adjoint system that is similar to the treatment described in [[5](/extras/refs#adjoint_cpr)].
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.CPRPreconditioner' href='#JutulDarcy.CPRPreconditioner'><span class="jlbinding">JutulDarcy.CPRPreconditioner</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CPRPreconditioner(p = default_psolve(), s = ILUZeroPreconditioner(); strategy = :quasi_impes, weight_scaling = :unit, update_frequency = 1, update_interval = :iteration, partial_update = true)
```


Construct a constrained pressure residual (CPR) preconditioner.

By default, this is a AMG-BILU(0) version (algebraic multigrid for pressure, block-ILU(0) for the global system).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/cpr.jl#L79-L86" target="_blank" rel="noreferrer">source</a></Badge>

</details>


The short version of the CPR preconditioner can be motivated by our test system:

$$r_n = \frac{\partial}{\partial t} ((1 - S_w) \rho_n \phi) + \nabla \cdot (\rho_n \vec{v}_n) - \rho_n q_n = 0,\\
r_w = \frac{\partial}{\partial t} (S_w \rho_w \phi) + \nabla \cdot (\rho_w \vec{v}_w) - \rho_w q_w = 0.$$

For simplicity, we assume that there is no gravity, source terms, or compressibility. Each equation can then be divided by their respective densities and summed up to produce a pressure equation:

$$r_p = \frac{\partial}{\partial t} ((1 - S_w) \phi) + \nabla \cdot \vec{v}_n + \frac{\partial}{\partial t} (S_w \phi) + \nabla \cdot  \vec{v}_w \\
= \frac{\partial}{\partial t} ((S_w - S_w) \phi) + \nabla \cdot (\vec{v}_n + \vec{v}_w) \\
= \nabla \cdot (\vec{v}_n + \vec{v}_w) \\
= - \nabla \mathbf{K}(k_{rw}/\mu_w + k_{rn}/\mu_n) \nabla p \\
= - \nabla \mathbf{K}\lambda_t \nabla p = 0$$

The final equation is the variable coefficient Poisson equation and is referred to as the incompressible pressure equation for a porous  media. We know that algebraic multigrid preconditioners (AMG) are highly efficient for linear systems made by discretizing this equation. The idea in CPR is to exploit this by constructing an approximate pressure equation that is suited for AMG inside the preconditioner.

Constructing the preconditioner is done in two stages:
1. First, weights for each equation is found locally in each cell that decouples the time derivative from the non-pressure variables. In the above example, this was the true IMPES weights (dividing by density). JutulDarcy supports analytical true IMPES weights for some systems, numerical true IMPES weights for all systems and quasi IMPES weights for all systems.
  
2. A pressure equation is formed by weighting each equation by the respective weights and summing. We then have two systems: The pressure system $r_p$ with scalar entries and the full system $r$ that has block structure.
  

During the linear solve, the preconditioner is then made up of two broad stages: First, a preconditioner is applied to the pressure part (typically AMG), then the full system is preconditioned (typically ILU(0)) after the residual has been corrected by the pressure estimate:
1. Form weighted pressure residual $r_p = \sum_i w_i r_i$.
  
2. Apply pressure preconditioer $M_p$: $\Delta p = M_p^{-1} r_p$.
  
3. Correct global residual $r^* = r - J P(\Delta p)$ where $P$ expands the pressure update to the full system vector, with zero entries outside the pressure indices.
  
4. Precondition the full system $\Delta x^* = M^{-1}r^*$
  
5. Correct the global update with the pressure to obtain the final update: $\Delta x = \Delta x^* + P(\Delta p)$
  

#### Multi model (porous medium with wells) {#Multi-model-porous-medium-with-wells}

If a model is a porous medium with wells, the same preconditioners can be used, but an additional step is required to incorporate the well system. In practical terms, this means that our linearized system is expanded to multiple linear systems:

$$J \Delta \mathbf{x} = \begin{bmatrix}
   J_{rr} & J_{rw} \\
   J_{wr} & J_{ww}
\end{bmatrix}
\begin{bmatrix}
\Delta \mathbf{x}_r \\
\Delta \mathbf{x}_w
\end{bmatrix}
 = 
\begin{bmatrix}
\mathbf{r}_r \\
\mathbf{r}_w
\end{bmatrix}$$

Here, $J_{rr}$ is the reservoir equations differentiated with respect to the reservoir primary variables, i.e. the Jacobian from the previous section. $J_{ww}$ is the well system differentiated with respect to the well primary variables. The cross terms, $J_{rw}$and $J_{wr}$, are the same equations differentiated with respect to the primary variables of the other system.

The well system is generally much smaller than the reservoir system and can be solved by a direct solver. We would like to reuse the block preconditioners defined for the base system. The approach we use is a [Schur complement](https://en.wikipedia.org/wiki/Schur_complement) approach to solve the full system. If we linearly eliminate the dependence of the reservoir equations on the well primary variables, we obtain the reduced system:

$$J \Delta \mathbf{x} = \begin{bmatrix}
   J_{rr} - J_{rw}J_{ww}^{-1}J_{wr} & 0 \\
   J_{wr} & J_{ww}
\end{bmatrix}
\begin{bmatrix}
\Delta \mathbf{x}_r \\
\Delta \mathbf{x}_w
\end{bmatrix}
 = 
\begin{bmatrix}
\mathbf{r}_r - J_{rw}J_{ww}^{-1}\mathbf{r}_w\\
\mathbf{r}_w
\end{bmatrix}$$

We can then solve the system in terms of the reservoir degrees of freedom where the system is a block linear system and we already have a working preconditioner:

$$\left(J_{rr} - J_{rw}J_{ww}^{-1}J_{wr}\right)\mathbf{x}_r = \mathbf{r}_r - J_{rw}J_{ww}^{-1}\mathbf{r}_w$$

Once that system is solved for $\mathbf{x}_r$, we can recover the well degrees of freedom $\mathbf{r}_w$ directly:

$$\mathbf{r}_w = J_{ww}^{-1}(\mathbf{r}_w - J_{wr}\mathbf{x}_r)$$

::: tip Efficiency of Schur complement

Explicitly forming the matrix $J_{rr} - J_{rw}J_{ww}^{-1}J_{wr}$ will generally lead to a lot of fill-in in the linear system. JutulDarcy instead uses the action of $J_{rr} - J_{rw}J_{ww}^{-1}J_{wr}$ as a linear operator from [LinearOperators.jl](https://github.com/JuliaSmoothOptimizers/LinearOperators.jl). This means that we must apply the inverse of the well system every time we need to compute the residual or action of the system matrix, but fortunately performing the action of the Schur complement is inexpensive as long as $J_{ww}$ is small and the factorization can be stored.

:::
