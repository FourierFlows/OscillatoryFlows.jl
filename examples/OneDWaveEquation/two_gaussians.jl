# # Two Gaussians
# 
# In this example, a Gaussian initial condition is let loose and two gaussians 
# that propagate in either direction emerge.

using
    OscillatoryFlows.OneDWaveEquation,
    Printf,
    Plots

# ## Build the problem

problem = Problem(; nx=256, Lx=2π, c=1, b=0, dt=0.01)

nothing # hide

# ## A Gaussian initial condition
#
# This initial condition has width $ \sqrt(2/50) = 1/5 $

ξ₀(x) = exp(-50x^2)

set_ξ!(problem, ξ₀)

# ## Run, and animate the results

anim = @animate 

plot(problem.grid.x, problem.vars.ξ,
      title = @sprintf("Two Gaussians, t = %.2f", problem.clock.t),
      ylims = (-0.5, 1.1), 
     legend = false))

for i = 1:50
    stepforward!(problem, 2)
    updatevars!(problem)

    plot(problem.grid.x, problem.vars.ξ,
          title = @sprintf("Two Gaussians, t = %.2f", problem.clock.t),
          ylims = (-0.5, 1.1), 
         legend = false))
end

mp4(anim, "two_gaussians.mp4", fps=8)
