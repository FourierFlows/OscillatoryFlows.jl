using OscillatoryFlows.OneDWaveEquation
using PyPlot

Lx = 2π
nx = 128
c = 1

problem = Problem(; nx=nx, Lx=Lx, c=c, dt=0.01)

ξ₀(x) = exp(-8x^2)

set_ξ!(problem, ξ₀)

close("all")
fig, axs = subplots()

updatevars!(problem)
plot(problem.grid.x, real.(problem.vars.ξ))

for i = 1:10
    stepforward!(problem, 1)
    updatevars!(problem)
    cla()
    plot(problem.grid.x, real.(problem.vars.ξ))
    pause(0.1)
end
