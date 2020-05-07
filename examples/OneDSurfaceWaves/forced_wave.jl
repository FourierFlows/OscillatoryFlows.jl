# # A forced Gaussian packet
# 
# In this examples, we illustrate the forcing of a Gaussian packet,
# by an idealized prescribed pressure forcing.

using 
    OscillatoryFlows.OneDSurfaceWaves,
    Printf,
    Statistics,
    PyPlot

# # Physical and numerical parameters

## Physical and numerical domain propeties
Lx = 2π     # Domain length
nx = 512    # Number of modes
 g = 1      # Gravitational acceleration
 h = 2      # Domain depth

## Wavenumber and implied frequency
 k = 32
 σ = √(g * k * tanh(k * h))
cg = g / 2σ * (k * h * sech(k * h)^2 + tanh(k * h))

## Modulation
 μ = 0.1
 ℓ = 1 / (μ * k)

## Forcing amplitude, duration, and beginning
ϵ₀ = 1e-3
a₀ = ϵ₀ / k
 T = 4 / σ
t₀ = 3T

## Initial condition
x₀ = 0.0 

## Safe time step
dt = 0.05 / σ

# ## Forcing definition

## Time-dependent amplitude
a(t) = a₀ / (g * t₀) * exp(-(t - t₀)^2 / (2 * T^2))

## Modulated resonant excitation
ϖ(x, t) = a(t) * exp(-(x - (x₀ + cg * t))^2 / (2 * ℓ^2)) * sin(k * x - σ * t)

grid = OneDGrid(nx, Lx, dealias=1/2)

problem = Problem(; grid=grid, dt=dt, g=g, h=h, ϖ=ϖ, order=4,
                    innerK=0.2, outerK=0.4)

# ## The velocity potential
potential = VelocityPotential(problem, nz=128)

## 2D meshes for plotting
X = repeat(potential.x, 1, potential.nz)
Z = repeat(potential.z, grid.nx, 1)

close("all")
fig, axs = subplots(nrows=2, sharex=true)

"""Plot the surface elevation and time."""
function makeplot!(problem)

    calculate!(potential)

    sca(axs[1]); cla()
    plot(k * problem.grid.x, problem.vars.s)
    title(@sprintf("\$ t / \\sigma = %.3f \$", problem.clock.t / σ))
    xlabel(L"k x"); ylabel(L"s")
    ylim(-5a₀, 5a₀)

    wmax = maximum(abs, potential.w)
    wlim = 2 * wmax / 100000

    sca(axs[2]); cla()
    pcolormesh(k * X, Z, potential.w, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

    pause(0.1)

    return nothing
end

## Plot initial condition
makeplot!(problem)

# # Run the problem
for i = 1:1000
    stepforward!(problem, 10)
    makeplot!(problem)
end
