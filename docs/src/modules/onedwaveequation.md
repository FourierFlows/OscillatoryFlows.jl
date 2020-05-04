# OneDWaveEquation Module

This module solves the one-dimensional wave equation:

```math
\xi_{tt} - c^2 \xi_{xx} = - \beta \xi_t
```

where ``\xi`` is the displacement associated with the wave field, ``c`` is
the wave speed, and ``\beta`` is a damping rate.

## Numerical formulation

Before solving the wave equation numerically, we 'diagonalize' the dispersion term by defining

```math
\chi = u + c \xi_x \\[1ex]
\psi = u - c \xi_x
```

To derive equations for ``\chi`` and ``\psi``, we start with a first-order formulation of the wave equation,

```math
\xi_t = u \\[1ex]
u_t = c^2 \xi_{xx} - \beta u
```

where ``u`` is the velocity, or the rate of change of the displacement.
We find

```math
\chi_t - c \chi_x = - \beta u \\[1ex]
\psi_t + c \psi_x = - \beta u
```

Taking the Fourier transform, we have
```math
\hat \chi = \hat u + \mathrm{i} c k \hat \xi \\[1ex]
\hat \psi = \hat u - \mathrm{i} c k \hat \xi 
```

We define

```math
\sigma \equiv c k
```

and note that

```math
  \hat u = \frac{1}{2} \left ( \hat \chi + \hat \psi \right ) \\[1ex]
\hat \xi = \frac{\mathrm{i}}{2 \sigma} \left ( \hat \psi - \hat \chi \right )
```

to obtain

```math
\hat \chi_t - \mathrm{i} \sigma \hat \chi = - \frac{\beta}{2} \left ( \hat \chi + \hat \psi \right ) \\[1ex]
\hat \psi_t + \mathrm{i} \sigma \hat \psi = - \frac{\beta}{2} \left ( \hat \chi + \hat \psi \right )
```

Using this form allows us to use the ETDRK4 time-stepper to solve the oscillatory part of this system exactly.
