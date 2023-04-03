# Solving the equations
By default, Jutul solves a system as a fully-coupled implicit system of equations.
## Newton's method
The standard way of solving a system of non-linear equations is by Newton's method (also known as Newton-Raphson's method). A quick recap: For a vector valued residual ``\mathbf{r}(x)`` of the primary variable vector ``\mathbf{x}`` we can defined a Newton update:
```math
\mathbf{x}^{k+1} = \mathbf{r}^k - J^{-1} \mathbf{r}(\mathbf{x}^k), \quad J_{ij} = \frac{\partial \mathbf{r}_i}{\partial \mathbf{x}_j}.
```
JutulDarcy solves systems that generally have both non-smooth behavior and physical constraints on the values for ``\textbf{x}``. For that reason, we modify Newton's method slightly:

```math
\mathbf{x}^{k+1} = \mathbf{r}^k + \omega(\Delta \mathbf{x})
```

Here, ``\oemga`` is a function that limits the variables so that they do not change too much (e.g. Appleyard chopping, limiting of pressure, saturation and composition updates) and that they are within the prescribed limits. There are also options for automated global dampening in the presence of convergence issues. The update is then defined from inverting the Jacobian:

```math
\Delta \mathbf{x} = -J^{-1} \mathbf{r}(\mathbf{x}^k), \quad J_{ij} = \frac{\partial \mathbf{r}_i}{\partial \mathbf{x}_j}.
```

Starting with ``\mathbf{x}^0``as some initial guess taken from the previous time-step, we can solve the system by iterating upon this loop.
## Linear solvers and linear systems
For most practical applications it is not feasible or efficient to invert the Jacobian. JutulDarcy uses preconditioned iterative solvers by default, but it is possible to use direct solvers as well when working with smaller models. The high level interface for setting up a reservoir model [`setup_reservoir_model`](@ref) has an optional `block_backend=true` keyword argument that determines the matrix format, and consequently the linear solver type to be used.

### Direct solvers
If `block_backend` is set to `false`, Jutul will assemble into the standard Julia CSC sparse matrix with `Float64` elements and Julia's default direct solver will be used. It is also possible to use other Julia solvers on this system, but the default preconditioners assume that block backend is enabled.

### Iterative solver 
If `block_backend` is set to `true`, Jutul will by default use a constrained-pressure residual (CPR) preconditioner for BiCGStab. Jutul relies on [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) for iterative solvers. 

#### Single model (only porous medium)
If the model is a single model (e.g. only a reservoir) the matrix format is a block-CSC matrix that combines Julia's builtin sparse matrix format with statically sized elements from the [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) package. If we consider the two-phase immiscible system from [Multi-phase, immiscible flow](@ref) we have a pair of equations ``R_n, R_w`` together with the corresponding primary variables pressure and first saturation ``p, S_n`` defined for all ``N_c`` cells. Let us simplify the notation a bit so that the subscripts of the primary variables are ``p, s`` and define a ``N_c \times N_c`` block Jacobian linear system where the entires are given by:
```math
J_{ij} = \begin{bmatrix}
   \left(\frac{\partial r_{n}}{\partial p}\right)_{ij} & \left(\frac{\partial r_{n}}{\partial s}\right)_{ij} \\
   \left(\frac{\partial r_{w}}{\partial p}\right)_{ij} & \left(\frac{\partial r_{w}}{\partial s}\right)_{ij} \end{bmatrix} = \begin{bmatrix}
   J_{np} & J_{ns} \\
   J_{wp} & J_{ws}
\end{bmatrix}_{ij}
```
This block system has several advantages:

 - We immediately get access to more powerful version of standard Julia preconditioners provided that all operations used are applicable for matrices and are applied in the right commutative order. For example, JutulDarcy uses the [ILUZero.jl](https://github.com/mcovalt/ILUZero.jl) package when a CSC linear system is preconditioned with incomplete LU factorization with zero fill-in.
 - Sparse matrix vector products are much more efficient as less indicies need to be looked up for each element wise multiplication.
 - Performing local reductions over variables is much easier when they are located in a local matrix.


##### CPR

#### Multi model (porous medium with wells)


### CSR backend and experimental features
`backend=:csr`
## References

CPR:

[Wallis, J.R. "Incomplete Gaussian Elimination as a Preconditioning for Generalized Conjugate Gradient Acceleration." Paper presented at the SPE Reservoir Simulation Symposium, San Francisco, California, November 1983](https://doi.org/10.2118/12265-MS)

[Cao, H., Tchelepi, H. A., Wallis, J., and H. Yardumian. "Parallel Scalable Unstructured CPR-Type Linear Solver for Reservoir Simulation." Paper presented at the SPE Annual Technical Conference and Exhibition, Dallas, Texas, October 2005](https://doi.org/10.2118/96809-MS)