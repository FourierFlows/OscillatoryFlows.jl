"""
    OneDSurfaceWaves

Solver for the potential flow equations in (x, z) beneath a 
free surface s(x, t). The system is one-dimensional because
the z-problem is solved analytically.
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
function Problem(T = Float64;
                       nx = 256,
                       Lx = 2π,
                        g = 9.81,
                        h = 1,
                        ϖ = zerofunction,
                     grid = OneDGrid(nx, Lx),
                    order = 3,
                   params = Params(T, g=g, h=h, ϖ=ϖ),
                  stepper = :RK4,
                       dt = 0.01,
                  #stepperkwargs...
                 )

      vars = Vars(order, grid)
       eqn = WaveEquation(grid)

    #return FourierFlows.Problem(eqn, stepper, dt, grid, vars, params, CPU();
    #                            stepperkwargs...)
    return FourierFlows.Problem(eqn, stepper, dt, grid, vars, params, CPU())
end

struct Params{T, P} <: AbstractParams
    g :: T 
    h :: T
    ϖ :: P
end

zerofunction(x, t) = 0

"""
    Params(T=Float64; g=9.81, h=1, ϖ=zerofunction)

Returns parameters for `OneDSurfaceWaves` problem with
`g`raviational acceleration domain `h`eight (or depth), and 
atmospheric pressure forcing, `ϖ` normalized by the density
of the fluid.
"""
Params(T=Float64; g=9.81, h=1, ϖ=zerofunction) = Params{T, typeof(ϖ)}(g, h, ϖ)

"""
    SurfaceWaveEquation(grid)

Returns an equation for surface waves on `grid`.
"""
SurfaceWaveEquation(grid) =
    FourierFlows.Equation(0, calcN!, grid, dims=(grid.nkr, 2))

"""
    PotentialPerturbationExpansion(T, order, grid)

Returns a Tuple of NamedTuples corresponding to 
grid values and Fourier coefficients of the velocity potential
evaluated at the surface on the one-dimensional `grid` in `x`.
The elements of the tuple are the terms in a perturbation expansion
of the velocity potential (and its derivatives) to `order`.
`T` is the floating point type of the arrays.

The members of the tuple are indexed by their order.
In other words, `ϕz[1]` and `ϕz[2]` are the O(ϵ¹) and O(ϵ²)
terms in the expansion of `ϕz`.

Each term in the expansion is a `NamedTuple`. The fields of the NamedTuple
are the vertical dervatives of the velocity potential.
For example,

    * `ϕ[1].o = ϕ[1][1]` is `ϕ₁`.
    * `ϕ[1].z = ϕ[1][2]` is the first z-derivative of `ϕ₁`
    * `ϕ[2].zzz` is the third z-derivative of `ϕ₂`

Each derivative contains both Fourier-space 'coefficients' and grid-space
arrays. Thus

    * `ϕ[3].zz.g` returns an array with grid space values for the
       second z-derivative of `ϕ₃`.

    * `ϕ[7].zzzz.c` returns a (complex) array with Fourier coefficients
       for the fourth z-derivative of `ϕ₇`.

There are `order` elements in the return tuple `ϕ`.
Each element `ϕ[m]` has `(order-m+1) elements.
Summary of the nested structure of `ϕ`:

    * `ϕ`: `Tuple` of `NamedTuple`s
    * `ϕ[m]`: `NamedTuple` of vertical derivatives
    * `ϕ[m].o`: `NamedTuple` of Fourier coefficients and grid-space values.
    * `ϕ[m].z`: `NamedTuple` of Fourier coefficients and grid-space values.
    * `ϕ[m].z.c`: Fourier coefficient array
    * `ϕ[m].z.g`: Grid space array

"""
function PotentialPerturbationExpansion(T, order, grid)

    ϕ̃ = []

    for m = 1:order
        derivative_names = [Symbol("z"^μ) for μ = 1:(order-m+1)]
        names = (:s, derivative_names...)

        fields = Tuple(
                       (c = zeros(Complex{T}, grid.nkr), g = zeros(T, grid.nx))
                       for μ = 1:(order-m+2)
                      )

        push!(ϕ̃, NamedTuple{names}(fields))
    end

    # Convert from abstractly-typed array to concretely-typed tuple
    return tuple(ϕ̃...)
end

"""
    inverse_transform!(ϕm, grid, N=1:length(ϕm))

Perform an inverse Fourier transform of the `N` fields
in `ϕm`, the O(ϵᵐ) term in a perturbation expansion of the
velocity potential.
"""
function inverse_transform!(ϕm, grid, N=1:length(ϕm))

    for n in N
        ∂ⁿz_ϕ = ϕm[n]
        ldiv!(∂ⁿz_ϕ.g, grid.rfftplan, ∂ⁿz_ϕ.c)
    end

    return nothing
end

"""
    forward_transform!(ϕm, grid, N=1:length(ϕm))

Perform a forward Fourier transform of the `N` fields
in `ϕm`, the O(ϵᵐ) term in a perturbation expansion of the
velocity potential.
"""
function forward_transform!(ϕm, grid, N=1:length(ϕm))

    for n in N
        ∂ⁿz_ϕ = ϕm[n]
        mul!(∂ⁿz_ϕ.c, grid.rfftplan, ∂ⁿz_ϕ.g)
    end

    return nothing
end

"""
    calculate_ϕm!(ϕ, m, s, grid) 

Computes the value of the mᵗʰ component of ϕz at the surface
using previously computed lower order vertical derivatives contained in ϕz,
according to the Taylor expansion developed in Dommermuth and Yue, 1987.

For `m = 2:4`, this does

    `ϕ₂ = - s * ϕ₁.z`    
    `ϕ₃ = - s * ϕ₂.z - s^2 / 2 * ϕ₁.zz`    
    `ϕ₄ = - s * ϕ₃.z - s^2 / 2 * ϕ₂.zz - s^3 / 6 * ϕ₁.zzz`    
"""
function calculate_ϕm!(ϕ, m, s, grid)

    @assert m > 1 "calculate_ϕm! is only correct for m > 1"

    ϕm = ϕ[m][1]

    @. ϕm.g = 0

    for n = 1:m-1
        @. ϕm.g -= s^n * ϕ[m - n][n + 1].g / factorial(n)
    end

    return nothing
end

"""
    ∂z_Ψ(m, h, k)

Compute the `m`ᵗʰ vertical derivative of the solution to
Laplace's equation in Fourier-space evaluated at `z=0`,
in a domain with depth `h` and with Fourier wavenumber `k`.

The mode is

    ``Ψ(z) = cosh(k * (z + h)) / cosh(k * h)``

The vertical derivatives at ``z=0`` are therefore

    *   `Ψz = k * tanh(k * h)`
    *  `Ψzz = k^2`
    * `Ψzzz = k^3 * tanh(k * h)`
"""
∂z_Ψ(m, h, k) = isodd(m) ? k^m * tanh(k * h) : k^m

"""
    compute_ϕz!(Σϕz, ϕzm, Φ̂, s, params, grid)

Computes the vertical gradient of the potential ϕ using a perturbation
expansion at the surface.
"""
function compute_ϕz!(Σϕz, ϕ, Φ̂, s, params, grid)
    @inbounds begin

        #####
        ##### First, we find ϕ₁ and its derivatives
        #####
        
        M = length(ϕ) # order of expansion
        ϕ₁ = ϕ[1] # Leading-order ϕ, ϕz, ϕzz, etc.

        # Obtain ϕ₁ and its derivatives in Fourier space
        for n = 0:M
            @. ϕ₁[n+1].c = Φ̂ * ∂z_Ψ(n, params.h, grid.kr)
        end

        # Transform ϕ₁ and its derivatives to grid space
        inverse_transform!(ϕ₁, grid)

        #####
        ##### Obtain ϕm ~ O(ϵᵐ) for m = 2, ..., M
        #####
        
        for m = 2:M
            # Set ϕ[m][1].g, the surface grid values of ϕ[m]
            calculate_ϕm!(ϕ, m, s, grid)

            # Convenient alias
            ϕm = ϕ[m]

            # Transform ϕm[1] to Fourier space
            mul!(ϕm[1].c, grid.rfftplan, ϕm[1].g)

            # Determine z-derivatives from ϕm[1].c
            for n = 2:length(ϕm)
                @. ϕm[n].c = ϕm[1].c * ∂z_Ψ(n, params.h, grid.kr)
            end

            # Transform the just-calculated z-derivative of ϕzm 
            # from coefficient space to grid space.
            inverse_transform!(ϕm, grid, 2:length(ϕm))
        end
    end

    # Given ϕ (and derivatives) at z=0, use a Taylor expansion 
    # to determine ϕz at z=s.
    surface_ϕz_from_expansion!(Σϕz, ϕ, s) 

    return nothing
end

"""
    surface_ϕz_from_expansion!(Σϕz, ϕ, s) 

Calculate the value of ϕz at z=s in terms of ϕ at z=0, and s, 
using a Taylor expansion.
"""
function surface_ϕz_from_expansion!(Σϕz, ϕ, s) 

    # Zero out Σϕz
    Σϕz .= 0

    # For each ϕm term...
    for m = 1:length(ϕ)
        ϕm = ϕ[m]

        # add every contribution to the Taylor series:
        for k = 0:length(ϕm) - 2
            ∂zᵏ⁺¹_ϕᵐ = ϕm[k + 2].g

            @. Σϕz += s^k / factorial(k) * ∂zᵏ⁺¹_ϕᵐ
        end
    end

    return nothing
end

struct Vars{P, A, B} <: AbstractVars
             ϕ :: P
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
function Vars(order, grid::AbstractGrid{T}) where {Dev, T}

    @devzeros CPU T grid.nx s ϖ Φ Φₓ sₓ ϕz aux_grid
    @devzeros CPU Complex{T} grid.nkr ϖ̂ ŝ aux_coeffs

    ϕ = PotentialPerturbationExpansion(T, order, grid)

    return Vars(ϕ, Φ, s, ϖ, Φₓ, sₓ, ϕz, aux_grid, ϖ̂, ŝ, aux_coeffs)
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

    ldiv!(vars.s, grid.rfftplan, ŝ)

    # The Dirichlet-to-Neumann map!
    compute_ϕz!(vars.ϕz, vars.ϕ, Φ̂, vars.s, params, grid)

    ######
    ###### Calculate Ns = ∂t s
    ######
    
    # Calculate (1 + sₓ²) ϕz
    ŝₓ = vars.aux_coeffs
    @views @. ŝₓ = im * grid.kr * ŝ
    ldiv!(vars.sₓ, grid.rfftplan, ŝₓ)

    @. vars.aux_grid = (1 + vars.sₓ^2) * vars.ϕz
    mul!(vars.aux_coeffs, grid.rfftplan, vars.aux_grid)

    @. Ns = vars.aux_coeffs

    # Calculate Φₓ sₓ
    Φ̂ₓ = vars.aux_coeffs
    @views @. Φ̂ₓ = im * grid.kr * Φ̂
    ldiv!(vars.Φₓ, grid.rfftplan, Φ̂ₓ)

    Φₓsₓ = vars.aux_grid
    Φₓsₓ_hat = vars.aux_coeffs
    @. Φₓsₓ = vars.Φₓ * vars.sₓ
    mul!(Φₓsₓ_hat, grid.rfftplan, Φₓsₓ)

    # Subtract the transform of Φₓsₓ
    @. Ns -= Φₓsₓ_hat

    ######
    ###### Calculate NΦ = ∂t Φ
    ######
    
    # Set NΦ = - g ŝ
    @views @. NΦ = - params.g * ŝ

    # Subtract forcing by ϖ
    @. vars.ϖ = params.ϖ(grid.x, t)
    mul!(vars.ϖ̂, grid.rfftplan, vars.ϖ)
    @views @. NΦ -= vars.ϖ̂

    # Calculate Φₓ²
    Φₓ² = vars.aux_grid
    Φₓ²_hat = vars.aux_coeffs
    @. Φₓ² = vars.Φₓ^2
    mul!(Φₓ²_hat, grid.rfftplan, Φₓ²)

    # Subtract transform of -Φₓ² / 2
    @. NΦ -= Φₓ²_hat / 2

    # Calculate (1 + sₓ²) ϕz²
    @. vars.aux_grid = (1 + vars.sₓ^2) * vars.ϕz^2
    mul!(vars.aux_coeffs, grid.rfftplan, vars.aux_grid)

    # Subtract (1 + sₓ²) ϕz² / 2
    @. NΦ += vars.aux_coeffs / 2
    
    @views dealias!(Ns, grid)
    @views dealias!(NΦ, grid)

    return nothing
end

"""
    updatevars!(vars, grid, sol)

Update `vars` on `grid` with the `sol`ution.
"""
function updatevars!(vars, params, grid, sol)
    ŝ = view(sol, :, 1)
    Φ̂ = view(sol, :, 2)

    compute_ϕz!(vars.ϕz, vars.ϕ, Φ̂, vars.s, params, grid)

    ldiv!(vars.s, grid.rfftplan, ŝ)
    ldiv!(vars.Φ, grid.rfftplan, Φ̂)

    @. vars.aux_coeffs = im * grid.kr * ŝ
    ldiv!(vars.sₓ, grid.rfftplan, vars.aux_coeffs)

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
    @. prob.vars.s = s
    @views mul!(prob.sol[:, 1], prob.grid.rfftplan, prob.vars.s)

    updatevars!(prob)

    return nothing
end

"""
    set_Φ!(prob, Φ)

Set the surface potential `Φ`.
"""
function set_Φ!(prob, Φ)
    @. prob.vars.Φ = Φ
    @views mul!(sol[:, 2], prob.grid.rfftplan, prob.vars.Φ)

    updatevars!(prob)

    return nothing
end

set_s!(prob, s::Function) = set_s!(prob, s.(prob.grid.x))
set_Φ!(prob, Φ::Function) = set_Φ!(prob, Φ.(prob.grid.x))

end # module
