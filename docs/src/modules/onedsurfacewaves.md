# OneDSurfaceWaves Module

This module solves the potential flow equations in ``x, z``.
The two prognostic variables are the surface elevation, ``s(x, t)``,
and the surface potential ``\Phi(x, t)`` defined via

```math
\Phi(x, t) \equiv phi(x, z=s(x, t), t) \, .
```

The velocity and pressure fields are related to ``\phi`` via

```math
\bm{u} \equiv \nabla \phi \, \qquad p \equiv - \phi_t \, .
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
\Phi_t = - g s - \frac{1}{2} \Phi^2_x - \frac{1}{2} \left ( 1 + s_x^2 \right ) \phi^2_z \, |_{z=s} - \varpi
```

where ``\varpi`` is the atmospheric pressure at ``z=s``.
