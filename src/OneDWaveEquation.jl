"""
    OneDWaveEquation

Solve the one dimensional wave equation with damping:

    ``∂t² ξ - c² ∂x² ξ = - β ∂t ξ``

To solve the problem, we 'diagonalize' dispersion by introducing

    ``χ = u + c ∂x ξ``
    ``ψ = u - c ∂x ξ``

In terms of these 'wave functions', the Fourier transform of the wave equation becomes

    ``∂t χ̂ - i σ χ̂ = - β/2 (ψ̂ + χ̂)``
    ``∂t ψ̂ + i σ ψ̂ = - β/2 (ψ̂ + χ̂)``

Note that 

    `` û = (ψ̂ + χ̂) / 2``
    `` ξ̂ = i (ψ̂ - χ̂) / 2 σ``
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
    T = typeof(grid.Lx)
    L = zeros(Complex{T}, grid.nkr, 2)

    # χ
    @. L[:, 1] =   im * σ(c, grid.kr)

    # ψ
    @. L[:, 2] = - im * σ(c, grid.kr)

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
    @devzeros Dev T grid.nx ξ u
    @devzeros Dev Complex{T} grid.nkr ξh uh
    return Vars(ξ, u, ξh, uh)
end

"""
    calcN!(N, sol, t, clock, vars, params, grid)

Calculate the nonlinear term for the 1D wave equation.
"""
function calcN!(N, sol, t, clock, vars, params, grid)
    @. N[:, 1] = - params.beta / 2 * (sol[:, 1] + sol[:, 2])
    @. N[:, 2] = - params.beta / 2 * (sol[:, 1] + sol[:, 2])
    return nothing
end

"""
    updatevars!(vars, grid, sol)

Update `vars` on `grid` with the `sol`ution.
"""
function updatevars!(vars, params, grid, sol)
    @. vars.uh = (sol[:, 2] + sol[:, 1]) / 2
    @. vars.ξh = (sol[:, 2] - sol[:, 1]) * im / (2 * σ(params.c, grid.kr))

    @inbounds vars.ξh[1] = 0

    ldiv!(vars.ξ, grid.rfftplan, vars.ξh)
    ldiv!(vars.u, grid.rfftplan, vars.uh)

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
    mul!(prob.vars.uh, prob.grid.rfftplan, prob.vars.u)

    # Set velocity u, and restore displacement ξ
    @. prob.sol[:, 1] = prob.vars.uh + im * σ(prob.params.c, prob.grid.kr) * prob.vars.ξh
    @. prob.sol[:, 2] = prob.vars.uh - im * σ(prob.params.c, prob.grid.kr) * prob.vars.ξh

    return nothing
end

"""
    set_ξ!(prob, ξ)

Set the wave displacement as the transform of `ξ`.
"""
function set_ξ!(prob, ξ)
    # Ensure that ξ̂ is saved:
    updatevars!(prob)

    # Set the velocity field
    prob.vars.ξ .= ξ
    mul!(prob.vars.ξh, prob.grid.rfftplan, prob.vars.ξ)

    # Set velocity u, and restore displacement ξ
    @. prob.sol[:, 1] = prob.vars.uh + im * σ(prob.params.c, prob.grid.kr) * prob.vars.ξh
    @. prob.sol[:, 2] = prob.vars.uh - im * σ(prob.params.c, prob.grid.kr) * prob.vars.ξh

    return nothing
end

set_u!(prob, u::Function) = set_u!(prob, u.(prob.grid.x))
set_ξ!(prob, ξ::Function) = set_ξ!(prob, ξ.(prob.grid.x))

end # module
