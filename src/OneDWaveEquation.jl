"""
    OneDWaveEquation

Solve the one dimensional wave equation with damping:

    ``∂t² ξ - c² ∂x² ξ = - β ∂t ξ``

To solve the problem, we use a complexified form of the Fourier transform of the
first-order formulation,

    ``∂t ξ = u``
    ``∂t u = c² ∂x² ξ - β u``

in which the Fourier-space wave function, ``φ̂ = σ ξ̂ + im û`` for ``σ = c k``, obeys

    ``∂t φ̂ + i σ φ̂ = - i β û``

Note that 

    `` û = im (  φ̂⋆ -  φ̂  ) / 2 ``
    `` ξ̂ =    (  φ̂  +  φ̂⋆ ) / 2σ ``
"""
module OneDWaveEquation

export
    Problem,
    updatevars!,
    set_u!,
    set_ξ!

using
    FFTW,
    Reexport

@reexport using FourierFlows

using LinearAlgebra: mul!, ldiv!

"""
    Problem(; parameters...)

Construct a wave propagation problem with steady or time-varying flow.
"""
function Problem(;
                       nx = 128,
                       Lx = 2π,
                     beta = 0,
                        c = 0.01,
                       dt = 0.01,
                  stepper = "ETDRK4",
                        T = Float64,
                 )

      grid = OneDGrid(CPU(), nx, Lx; T=T)
    params = Params{T}(c, beta)
      vars = Vars(CPU(), grid)
       eqn = WaveEquation(c, grid)

    return FourierFlows.Problem(eqn, stepper, dt, grid, vars, params, CPU())
end

struct Params{T} <: AbstractParams
       c :: T
    beta :: T
end

# The dispersion relation
σ(c, k) = c * k

"""
    WaveEquation(beta, grid)

Returns a wave equation with damping beta, on grid.
"""
function WaveEquation(c, grid)
    L = @. im * σ(c, grid.k)
    return FourierFlows.Equation(L, calcN!, grid)
end

struct Vars{A, B} <: AbstractVars
     ξ :: A
     u :: A
    ξh :: B
    uh :: B
end

"""
    Vars(dev, grid)

Returns the vars for a one-dimensional wave equation problem problem on `grid`.
"""
function Vars(::Dev, grid::AbstractGrid{T}) where {Dev, T}
    @devzeros Dev Complex{T} grid.nx ξ u
    @devzeros Dev Complex{T} grid.nk ξh uh
    return Vars(ξ, u, ξh, uh)
end

"""
    calcN!(N, sol, t, clock, vars, params, grid)

Calculate the nonlinear term for the 1D wave equation.
"""
function calcN!(N, sol, t, clock, vars, params, grid)
    @. N = params.beta / 2 * (conj(sol) - sol)
    return nothing
end

"""
    updatevars!(vars, grid, sol)

Update `vars` on `grid` with the `sol`ution.
"""
function updatevars!(vars, params, grid, sol)
    @. vars.ξh =      ( sol       + conj(sol) ) / (2 * σ(grid.k, params.c))
    @. vars.uh = im * ( conj(sol) - sol       ) / 2

    @inbounds vars.ξh[1] = 0

    ldiv!(vars.ξ, grid.fftplan, vars.ξh)
    ldiv!(vars.u, grid.fftplan, vars.uh)

    return nothing
end

"""
    updatevars!(prob)

Update `prob.vars` in the one-dimensional wave `prob`lem
"""
updatevars!(prob) = updatevars!(prob.vars, prob.params, prob.grid, prob.sol)

"""
    set_u!(prob, u)

Set the wave velocity as the transform of `u`.
"""
function set_u!(prob, u)
    # Ensure that ξ̂ is saved:
    updatevars!(prob)

    # Set the velocity field
    prob.vars.u .= u
    mul!(prob.vars.uh, prob.grid.fftplan, prob.vars.u)

    # Set velocity u, and restore displacement ξ
    @. prob.sol = im * prob.vars.uh
    @. prob.sol += σ(prob.params.c, prob.grid.k) * prob.vars.ξh

    return nothing
end

"""
    set_ξ!(prob, ξ)

Set the wave displacement as the transform of `ξ`.
"""
function set_ξ!(prob, ξ)
    # Ensure that û is saved:
    updatevars!(prob)

    # Set the displacement field
    prob.vars.ξ .= ξ
    mul!(prob.vars.ξh, prob.grid.fftplan, prob.vars.ξ)

    # Set displacement and restore velocity
    @. prob.sol = im * prob.vars.uh
    @. prob.sol += σ(prob.params.c, prob.grid.k) * prob.vars.ξh

    return nothing
end

set_u!(prob, u::Function) = set_u!(prob, u.(prob.grid.x))
set_ξ!(prob, ξ::Function) = set_ξ!(prob, ξ.(prob.grid.x))

end # module
