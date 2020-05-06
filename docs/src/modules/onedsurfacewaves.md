# OneDSurfaceWaves

This module solves the potential flow equations in ``x, z``
beneath a free surface at ``z = s(x, t)``.
The numerical method uses a Taylor expansion of the boundary condition 
at the free surface due to 
[Dommermuth and Yue (1987)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/highorder-spectral-method-for-the-study-of-nonlinear-gravity-waves/89E8FAC1EA7D3FB5F70D37ED4A2B2054).

## Governing equations

The two prognostic variables are the surface elevation, ``s(x, t)``,
and the surface potential ``\Phi(x, t)`` defined via

```math
\Phi(x, t) \equiv \phi(x, z=s(x, t), t) \, .
```

The velocity and pressure fields are related to ``\phi`` via

```math
\bm{u} \equiv \nabla \phi \, \qquad p \equiv - \phi_t - \frac{1}{2} \left | \bnabla \phi \right |^2 \, .
```

The surface elevation obeys

```math
s_t = \left ( 1 + s_x^2 \right ) \phi_z \, |_{z=s} - \Phi_x s_x \, .
```

where ``\phi_z \, |_{z=s}`` is the vertical gradient of the velocity potential --- the vertical velocity --
at ``z=s``.

The dynamic boundary condition on pressure at the surface of the ocean yields
an evolution equation for ``\Phi``:

```math
\Phi_t = - g s - \frac{1}{2} \Phi^2_x + \frac{1}{2} \left ( 1 + s_x^2 \right ) \phi^2_z \, |_{z=s} - \varpi
```

where ``\varpi`` is the atmospheric pressure at ``z=s``.

## Algorithm to obtain ``\phi_z \, |_{z=s}``

The main difficulty is finding the vertical velocity at the surface, ``\phi_z \, |_{z=s}``.

For this we introduce a perturbation expansion of ``\phi`` at ``z=s``,

```math
\phi(x, z, t) = \sum_{m=1}^M \phi_m(x, z, t) \, ,
```

where each ``\phi_m`` is small than ``\phi_{m-1}`` by ``O(\epsilon)``, where
``\epsilon`` is the surface slope (which must be small for the validity of this algorithm).

```math
\begin{align}
\Phi(x, t) &= \phi \at{z=s} \\
           & \approx \big [   \overbrace{\overstrut{2ex} \phi_1 }^{O(\ep)}
                            + \overbrace{\overstrut{2ex} \phi_2 + s \phi_{1z}}^{O \left ( \ep^2 \right )}
                            + \overbrace{\overstrut{2ex} \phi_3 + s \phi_{2z} + \frac{1}{2} s^2 \phi_{1zz}}^{O \left ( \ep^3 \right )} + \cdots \, \big ]_{z=0} \\
           & = \sum_{m=1}^M \sum_{n=0}^{M-m} \frac{s^n}{n!} \partial_z^n \phi_m \at{z=0} + O \left ( \epsilon^M \right ) \per
\end{align}
```


