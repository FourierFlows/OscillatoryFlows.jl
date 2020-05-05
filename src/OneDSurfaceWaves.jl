"""
    OneDSurfaceWaves
"""
module OneDSurfaceWaves

export
    Problem,
    updatevars!,
    set_s!,
    set_Φ!

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
                    order = 3,
                        g = 9.81,
                        h = 1,
                       dt = 0.01,
                  stepper = "RK4",
                        T = Float64,
                 )

      grid = OneDGrid(CPU(), nx, Lx; T=T)
    params = Params(T, g=g, h=h)
      vars = Vars(CPU(), grid, order)
       eqn = WaveEquation(grid)

    return FourierFlows.Problem(eqn, stepper, dt, grid, vars, params, CPU())
end

struct Params{T, P} <: AbstractParams
    g :: T 
    h :: T
    ϖ :: P
end

zerofunction(x, t) = 0

Params(T=Float64; g=9.81, h=1, ϖ=zerofunction) = Params{T}(g, h, ϖ)

"""
    WaveEquation(beta, grid)

Returns a wave equation with damping beta, on grid.
"""
function WaveEquation(grid)
    T = typeof(grid.Lx)
    L = zeros(Complex{T}, grid.nkr, 2)

    return FourierFlows.Equation(L, calcN!, grid)
end

function PotentialPerturbationExpansion(T, order)

    ϕ̃z = []

    for m = 1:order
        names = Tuple(Symbol("z"^μ) for μ = 1:(order-m))

        fields = Tuple(
                       (c = zeros(Complex{T}, grid.nkr), g = zeros(T, grid.nx))
                       for μ = 1:(order-m)
                      )

        push!(ϕ̃z, NamedTuple{names}(fields))
    end

    return Tuple(ϕ̃z...)
end

function inverse_transform!(ϕzm, grid, N=1:length(ϕzm))

    for n in N
        ∂ⁿz_ϕ = ϕzm[n]
        ldiv!(∂ⁿz_ϕ.g, grid.rfftplan, ∂ⁿz_ϕ.c)
    end

    return nothing
end

function forward_transform!(ϕzm, grid, N=1:length(ϕzm))

    for n in N
        ∂ⁿz_ϕ = ϕzm[n]
        mul!(∂ⁿz_ϕ.c, grid.rfftplan, ∂ⁿz_ϕ.g)
    end

    return nothing
end

∂z_Ψ(m, h, k) = isodd(m) ? k^m * tanh(k * h) : k^m

"""
    set_surface_ϕz!(ϕzm, m, ϕz, s, grid) 

Computes the value of the mᵗʰ component of ϕz at the surface
using previously computed lower order vertical derivatives contained in ϕz,
according to the Taylor expansion developed in Dommermuth and Yue, 1987.
"""
function set_surface_ϕz!(ϕzm, m, ϕz, s, grid)
    @. ϕzm = 0

    for n = (m-1) : -1 : 1
        @. ϕzm.g += s^n * ϕz[n] / factorial(n)
    end

    mul!(ϕzm.c, grid.rfftplan, ϕzm.g)

    return nothing
end

"""
    compute_ϕz!(Σϕz, ϕzm, Φ̂, s, params, grid)

Computes the vertical gradient of the potential ϕ using a perturbation
expansion at the surface.
"""
function compute_ϕz!(Σϕz, ϕzm, Φ̂, s, params, grid)
    M = length(ϕz)
    ϕz1 = ϕz[1]

    for n = 1:M
        @. ϕz1[n].c = Φ̂ * ∂z_Ψ(n, params.h, grid.kr)
    end

    inverse_transform!(ϕz1, grid)

    for m = 2:M
        ϕzm = ϕz[m]
        set_surface_ϕz!(ϕzm, m, ϕz, s, grid)

        for n = 2:(M-m)
            @. ϕzm[n].c = ϕzm[1].c / ∂z_Ψ(1, params.h, grid.kr) * ∂z_Ψ(m, params.h, grid.kr)
        end

        inverse_transform!(ϕzm, grid, 2:(M-m))
    end

    Σϕz .= 0
    for m = 1:M
        @. Σϕz += ϕz[1]
    end

    return nothing
end

struct Vars{P, A, B} <: AbstractVars
           ϕzm :: P
             Φ :: A
             s :: A
             ϖ :: A
            Φₓ :: A
            sₓ :: A
            ϕz :: A
      aux_grid :: A
             ϖ̂ :: B
             ŝ :: B
    aux_coeffs :: B
end

"""
    Vars(dev, grid)

Returns the vars for a one-dimensional free surface problem on `grid`.
"""
function Vars(::Dev, order, grid::AbstractGrid{T}) where {Dev, T}

    @devzeros Dev T grid.nx s ϖ Φ Φₓ sₓ ϕz aux_grid
    @devzeros Dev Complex{T} grid.nkr ϖ̂ ŝ aux_coeffs

    ϕzm = PotentialPerturbationExpansion(T, order)

    return Vars(ϕzm, Φ, s, ϖ, Φₓ, sₓ, ϕz, aux_grid, ϖ̂, ŝ, aux_coeffs)
end

"""
    calcN!(N, sol, t, clock, vars, params, grid)

Calculate the nonlinear term for the 1D wave equation.
"""
function calcN!(N, sol, t, clock, vars, params, grid)
    ŝ = view(sol, :, 1)
    Φ̂ = view(sol, :, 2)

    Ns = view(N, :, 1)
    NΦ = view(N, :, 2)

    ldiv!(vars.Φ, grid.rfftplan, ŝ)
    ldiv!(vars.s, grid.rfftplan, Φ̂)

    compute_ϕz!(vars.ϕz, vars.ϕzm, ŝ, vars.s, params, grid)

    ######
    ###### Calculate RHS for s
    ######
    
    # Calculate Φₓ sₓ
    Φ̂ₓ = vars.aux_coeffs
    @views @. Φ̂ₓ = im * grid.kr * ŝ
    ldiv!(vars.Φₓ, grid.rfftplan, Φ̂ₓ)

    ŝₓ = vars.aux_coeffs
    @views @. ŝₓ = im * grid.kr * Φ̂
    ldiv!(vars.sₓ, grid.rfftplan, ŝₓ)

    @. vars.aux_grid = vars.Φₓ * vars.sₓ
    mul!(vars.aux_coeffs, grid.rfftplan, vars.aux_grid)

    @. Ns = - vars.aux_coeffs

    # Calculate (1 + sₓ²) ϕz
    @. vars.aux_grid = (1 + vars.sₓ^2) * ϕz
    mul!(vars.aux_coeffs, grid.rfftplan, vars.aux_grid)

    @. Ns += vars.aux_coeffs

    ######
    ###### Calculate RHS for Φ
    ######
    
    # Calculate (1 + sₓ²) ϕz²
    @. vars.aux_grid = (1 + vars.sₓ^2) * ϕz^2
    mul!(vars.aux_coeffs, grid.rfftplan, vars.aux_grid)

    # Add -(1 + sₓ²) ϕz²
    @. NΦ = vars.aux_coeffs / 2

    # Add - g ŝ
    @. NΦ = - params.g * ŝ

    # Calculate Φₓ²
    @. vars.aux_grid = vars.Φₓ^2
    mul!(vars.aux_coeffs, grid.rfftplan, vars.aux_grid)

    # Add transform of -Φₓ² / 2
    @. NΦ = - vars.aux_coeffs / 2

    # Add forcing by ϖ
    vars.ϖ .= params.ϖ.(grid.x, t)
    mul!(vars.ϖ̂, grid.rfftplan, vars.ϖ)
    @. NΦ = - vars.ϖ̂

    return nothing
end

"""
    updatevars!(vars, grid, sol)

Update `vars` on `grid` with the `sol`ution.
"""
function updatevars!(vars, params, grid, sol)
    ŝ = view(sol, :, 1)
    Φ̂ = view(sol, :, 2)

    compute_ϕz!(vars.ϕz, vars.ϕzm, ŝ, vars.s, params, grid)

    @views ldiv!(vars.Φ, grid.rfftplan, ŝ)
    @views ldiv!(vars.s, grid.rfftplan, Φ̂)

    return nothing
end

"""
    updatevars!(prob)

Update `prob.vars` in the one-dimensional wave `prob`lem
"""
updatevars!(prob) = updatevars!(prob.vars, prob.params, prob.grid, prob.sol)

"""
    set_s!(prob, u)

Set the surface displacement `s`.
"""
function set_s!(prob, s)
    ŝ = view(sol, :, 1)
    Φ̂ = view(sol, :, 2)

    prob.vars.s .= s
    @views mul!(ŝ, prob.grid.rfftplan, prob.vars.s)

    updatevars!(prob)

    return nothing
end

"""
    set_Φ!(prob, Φ)

Set the surface potential `Φ`.
"""
function set_Φ!(prob, Φ)
    ŝ = view(sol, :, 1)
    Φ̂ = view(sol, :, 2)

    prob.vars.Φ .= Φ
    @views mul!(Φ̂, prob.grid.rfftplan, prob.vars.Φ)

    return nothing
end

set_s!(prob, u::Function) = set_s!(prob, s.(prob.grid.x))
set_Φ!(prob, ξ::Function) = set_Φ!(prob, Φ.(prob.grid.x))

end # module
