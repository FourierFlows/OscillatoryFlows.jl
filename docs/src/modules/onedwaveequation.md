# OneDWaveEquation Module

This module solves the one-dimensional wave equation:

```math
\xi_{tt} - c^2 \xi_{xx} = - \beta \xi_t
```

where ``\xi`` is the displacement associated with the wave field, ``c`` is
the wave speed, and ``\beta`` is a damping rate.

## Numerical formulation

To solve the wave equation, we use a complexification of the first-order formulation of the 
wave equation,

```math
\xi_t = u \\[1ex]
u_t = c^2 \xi_{xx} - \beta u
```

where ``u`` is the velocity, or the rate of change of the displacement.
Taking the Fourier transform of ``\xi`` and ``u`` and forming ``\sigma \hat \xi_t + \mathrm{i} \hat u_t``, 
where ``\sigma \equiv c k`` for each Fourier wavenumber, leads to

```math
\partial_t \left ( \sigma \hat \xi + \mathrm{i} \hat u \right ) = 
    - \mathrm{i} \sigma \left ( \sigma \hat \xi + \mathrm{i} \hat u \right ) - \beta \hat u \, .
```

We thus define the Fourier space "wave function"

```math
\hat \varphi = \sigma \hat \xi + \mathrm{i} u \, ,
```

which obeys

```math
\partial_t \hat \varphi + \mathrm{i} \sigma \hat \varphi = \frac{\beta}{2} \left ( \hat \varphi^\star - \hat \varphi \right ) \, .
```

Using this form allows us to use the ETDRK4 time-stepper to solve the oscillatory part of this system exactly.

## Recovering the primitive variables

The primitive variables ``\xi`` and ``u`` are recovered via

```math
\hat \xi = \frac{1}{2 \sigma}  \left ( \hat \varphi       + \hat \varphi^\star \right ) \\[2ex]
\hat u = \frac{\mathrm{i}}{2} \left ( \hat \varphi^\star - \hat \varphi       \right )
```
The resulting equation for ``\varphi`` is
