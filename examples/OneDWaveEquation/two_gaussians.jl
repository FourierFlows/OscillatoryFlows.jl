# # Two Gaussians
# 
# In this example, a Gaussian initial condition is let loose and two gaussians 
# that propagate in either direction emerge.

using
    OscillatoryFlows.OneDWaveEquation,
    Printf,
    Plots

# ## Build the problem

problem = Problem(; nx=256, Lx=2π, c=1, β=0, dt=0.01)

nothing # hide

# ## A Gaussian initial condition
#
# This initial condition is a Gaussian with variance $1/10$.

σ = 1/10 # variance
ξ₀(x) = exp( - x^2 / (2σ^2) )

set_ξ!(problem, ξ₀)

function makeplot(problem)
    p = plot(problem.grid.x, problem.vars.ξ,
          ylims = (-0.5, 1.1),
          title = @sprintf("Two Gaussians, t = %.2f", problem.clock.t),
         xlabel = "x",
         ylabel = "ξ",
         legend = false)
    return p
end

nothing # hide

# ## Run, and animate the results

anim = @animate for i = 1:51
    makeplot(problem)   # plot before stepforward!() to get frame with initial condition
    stepforward!(problem, 2)
    updatevars!(problem)
end

mp4(anim, "two_gaussians.mp4", fps=12)
