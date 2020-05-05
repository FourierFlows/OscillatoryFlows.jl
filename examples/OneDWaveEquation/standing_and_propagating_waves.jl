# # Standing and propagating waves
# 
# In this example, we illustrate the motion of both standing and propagating
# waves using `OscillatoryFlows`' `OneDWaveEquation` module.

using 
    OscillatoryFlows.OneDWaveEquation,
    Printf,
    PyPlot

# ## Parameters

Lx = 2π     # domain length
nx = 256    # number of grid points
 c = 1      # wave phase speed
dt = 0.01   # time-step

# ## Setting up the problems
   standing_wave_problem = Problem(; nx=nx, Lx=Lx, c=c, beta=0, dt=dt)
propagating_wave_problem = Problem(; nx=nx, Lx=Lx, c=c, beta=0, dt=dt)

# ## Initial conditions
# ξ(x) = cos(k * x - k * c * t)
# u(x) = k * c * sin(k * x - k * c * t)
k = 8
ξ₀(x) = cos(k * x)
u₀(x) = k * c * sin(k * x)

# ### Standing wave initial condition

set_ξ!(standing_wave_problem, ξ₀)

# ### Propagating wave initial condition

set_ξ!(propagating_wave_problem, ξ₀)
set_u!(propagating_wave_problem, u₀)


# ## Visualization

close("all")
fig, axs = subplots(nrows=2, sharey=true, sharex=true)

function makeplot!(axs, standing_wave_problem, propagating_wave_problem)

    updatevars!(standing_wave_problem)
    updatevars!(propagating_wave_problem)
    
    sca(axs[1]); cla()
    plot(standing_wave_problem.grid.x, real.(standing_wave_problem.vars.ξ))
    title(@sprintf("Standing wave, \$ t = %.2f \$", standing_wave_problem.clock.t))

    sca(axs[2]); cla()
    plot(propagating_wave_problem.grid.x, real.(propagating_wave_problem.vars.ξ))
    title(@sprintf("Propagating wave, \$ t = %.2f \$", propagating_wave_problem.clock.t))

    ylim(-1.1, 1.1)
    tight_layout()
    pause(0.05)

    return nothing
end

makeplot!(axs, standing_wave_problem, propagating_wave_problem)

# ## Time-step

for i = 1:100
    stepforward!(standing_wave_problem, 1)
    stepforward!(propagating_wave_problem, 1)

    makeplot!(axs, standing_wave_problem, propagating_wave_problem)
end
