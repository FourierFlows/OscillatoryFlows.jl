# # Two Gaussians
# 
# In this example, a gaussian initial condition is let loose and two gaussians 
# that propagate in either direction emerge.

using
    OscillatoryFlows.OneDWaveEquation,
    Printf,
    PyPlot

Lx = 2π
nx = 256
c = 1

problem = Problem(; nx=nx, Lx=Lx, c=c, beta=0, dt=0.01)

ξ₀(x) = exp(-64x^2)

set_ξ!(problem, ξ₀)

close("all")
fig, axs = subplots()

updatevars!(problem)
plot(problem.grid.x, real.(problem.vars.ξ))

for i = 1:100
    stepforward!(problem, 1)
    updatevars!(problem)
    cla()
    plot(problem.grid.x, real.(problem.vars.ξ))
    title(@sprintf("Two Gaussians, \$ t = %.2f \$", problem.clock.t))
    ylim(-0.5, 1.1)
    pause(0.1)
end
