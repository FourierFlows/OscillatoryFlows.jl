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
nx = 128
dt = 0.01

   linear_problem = Problem(; nx=nx, Lx=Lx, dt=dt, g=1, h=1, order=1,
                           stepper=:FilteredRK4)
nonlinear_problem = Problem(; nx=nx, Lx=Lx, dt=dt, g=1, h=1, order=3,
                           stepper=:FilteredRK4)

ℓ = 0.05

nonlinear_a = 1e-3
   linear_a = 1e-5
nonlinear_s₀(x) = nonlinear_a * exp(-x^2 / (2ℓ^2))
   linear_s₀(x) = linear_a * exp(-x^2 / (2ℓ^2))

# # Gaussian initial condition

set_s!(linear_problem,       linear_s₀)
set_s!(nonlinear_problem, nonlinear_s₀)

close("all")
fig, axs = subplots(nrows=2)

function makeplot!(linear_problem, nonlinear_problem)

    updatevars!(linear_problem)
    updatevars!(nonlinear_problem)

    sca(axs[1]); cla()
    plot(linear_problem.grid.x, linear_problem.vars.s / linear_a)
    plot(nonlinear_problem.grid.x, nonlinear_problem.vars.s / nonlinear_a)
    title(@sprintf("\$ t = %.1f \$", linear_problem.clock.t))
    ylim(-1.1, 1.1)

    sca(axs[2]); cla()
    plot(linear_problem.grid.x, linear_problem.vars.sₓ)
    plot(nonlinear_problem.grid.x, nonlinear_problem.vars.sₓ)

    pause(0.1)

    return nothing
end

makeplot!(linear_problem, nonlinear_problem)

for i = 1:100
    stepforward!(linear_problem, 10)
    stepforward!(nonlinear_problem, 10)
    makeplot!(linear_problem, nonlinear_problem)
end
