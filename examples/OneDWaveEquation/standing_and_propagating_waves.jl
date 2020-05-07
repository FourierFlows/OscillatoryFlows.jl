# # Standing and propagating waves
# 
# In this example, we illustrate the motion of both standing and propagating
# waves using `OscillatoryFlows`' `OneDWaveEquation` module.

using 
    OscillatoryFlows.OneDWaveEquation,
    Printf,
    Plots

# ## Parameters
#
# Here, we choose a domain size, numerical resolution,
# wave phase speed, and the time-step.

Lx = 2π     # domain length
nx = 256    # number of grid points
 c = 1      # wave phase speed
dt = 0.01   # time-step

nothing # hide

# ## Set up two problems
#
# We set up two problems: one to simulate a standing wave, and
# another to simulate a propagating wave.

   standing_problem = Problem(; nx=nx, Lx=Lx, c=c, beta=0, dt=dt)
propagating_problem = Problem(; nx=nx, Lx=Lx, c=c, beta=0, dt=dt)

nothing # hide

# ## Initial conditions
#
# For our initial conditions we use a sinusoid with wavenumber=8.
#
# The initial velocity for the standing wave is zero.
#
# We determine the initial velocity for the propagating wave 
# using the dispersion relation. For t > 0, the propagating
# wave solution has the displacement
#
# $ ξ(x, t) = \cos(k x - k c t) $
#
# This implies that the velocity of the standing wave,
# $u = ∂_t ξ$, is
#
# $ u(x, t) = k * c * \sin(k * x - k * c * t) $
#
# Taking $t=0$ determines the initial conditions, 
# $ξ(x, t=0)$ and $u(x, t=0)$.

## Wavenumber
k = 8

## Functions that define the initial conditions
ξ₀(x) = cos(k * x)
u₀(x) = k * c * sin(k * x)

## Set standing wave initial condition
set_ξ!(standing_problem, ξ₀)

## Set propagating wave initial condition
set_ξ!(propagating_problem, ξ₀)
set_u!(propagating_problem, u₀)

# ## Time-step and animate
#
# We're finally ready to time-step our problem forward.
# Along the way, we create an animation to visualize the solution.

anim = @animate for i = 1:25
    stepforward!(standing_problem, 4)
    stepforward!(propagating_problem, 4)

    updatevars!(standing_problem)
    updatevars!(propagating_problem)

    standing_plot = 
        plot(standing_problem.grid.x, standing_problem.vars.ξ,
             title = @sprintf("Standing wave, t = %.2f", standing_problem.clock.t),
             xlabel = "x", ylabel = "s", ylims = (-1, 1))
             
    propagating_plot = 
        plot(propagating_problem.grid.x, propagating_problem.vars.ξ,
             title = @sprintf("Propagating wave, t = %.2f ", standing_problem.clock.t),
             xlabel = "x", ylabel = "s", ylims = (-1, 1))

    plot(standing_plot, propagating_plot, layout = (2, 1), legend = false)
end

mp4(anim, "standing_propagating_waves.mp4", fps=8) # hide
