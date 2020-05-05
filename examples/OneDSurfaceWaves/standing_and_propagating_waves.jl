# # Standing and propagating waves
# 
# In this examples, we illustrate the motion of both standing
# and propagating waves using `OscillatoryFlows`' `OneDWaveEquation` module.

using 
    OscillatoryFlows.OneDSurfaceWaves,
    Printf,
    PyPlot

Lx = 2π
nx = 256
dt = 0.01

problem = Problem(; nx=nx, Lx=Lx, dt=dt)

a = 0.01
k = 8
s₀(x) = a * cos(k * x)

# # Standing wave initial condition

set_s!(problem, s₀)

close("all")
fig, axs = subplots()

function makeplot!(problem)

    cla()

    updatevars!(problem)
    
    plot(problem.grid.x, problem.vars.s)
    pause(0.05)

    return nothing
end

makeplot!(problem)

for i = 1:100
    stepforward!(problem, 1)
    makeplot!(problem)
end
