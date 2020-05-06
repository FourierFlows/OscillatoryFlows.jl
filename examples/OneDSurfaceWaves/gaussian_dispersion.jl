# # Gaussian dispersion
# 
# In this examples, we illustrate how a Gaussian initial surface elevation
# disperses into component waves that propagate to ±x, using 
# the OscillatoryFlows.OneDSurfaceWaves module.

using 
    OscillatoryFlows.OneDSurfaceWaves,
    Printf,
    PyPlot

Lx = 2π
nx = 1024
dt = 0.1

problem = Problem(; nx=nx, Lx=Lx, dt=dt, g=1, h=1, order=1)

a = 0.01
ℓ = 0.01
s₀(x) = a * exp(-x^2 / (2ℓ^2))

# # Gaussian initial condition

set_s!(problem, s₀)

close("all")
fig, axs = subplots()

function makeplot!(problem)

    updatevars!(problem)

    cla()
    plot(problem.grid.x, problem.vars.s)
    title(@sprintf("\$ t = %.1f \$", problem.clock.t))
    ylim(-1.1a, 1.1a)
    pause(0.1)

    return nothing
end

makeplot!(problem)

for i = 1:100
    stepforward!(problem, 1)
    makeplot!(problem)
end
