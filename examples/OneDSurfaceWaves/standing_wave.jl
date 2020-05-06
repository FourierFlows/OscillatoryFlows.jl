# # Standing and propagating waves
# 
# In this examples, we illustrate the motion of both standing
# and propagating waves using `OscillatoryFlows`' `OneDWaveEquation` module.

using 
    OscillatoryFlows.OneDSurfaceWaves,
    Printf,
    PyPlot

Lx = 2π
nx = 64 
dt = 0.1

problem = Problem(; nx=nx, Lx=Lx, dt=dt, g=1, h=1, order=1)

a = 0.01
k = 1
s₀(x) = a * cos(k * x)

# # Standing wave initial condition

set_s!(problem, s₀)

close("all")
fig, axs = subplots()

function makeplot!(problem)

    updatevars!(problem)

    cla()
    plot(problem.grid.x, problem.vars.s)
    title("$(problem.clock.t)")
    ylim(-1.1a, 1.1a)
    pause(0.1)

    return nothing
end

makeplot!(problem)

for i = 1:100
    stepforward!(problem, 1)
    makeplot!(problem)
end
