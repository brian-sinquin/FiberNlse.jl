# Theoretical Background

The propagation of ultrashort pulses in nonlinear waveguides is described by the Generalized Nonlinear Schrödinger Equation (GNLSE). This section presents both the fundamental theory and its practical implementation in FiberNlse.jl.

## The GNLSE

The implemented form of the GNLSE is:

```math
\frac{\partial A}{\partial z} = -\frac{\alpha}{2}A + i\sum_{k\geq 2}\frac{\beta_k}{k!}\left(i\frac{\partial}{\partial t}\right)^k A + i\gamma\left(1 + \frac{i}{\omega_0}\frac{\partial}{\partial t}\right)\left[|A|^2A(1-f_R) + A\int_{-\infty}^{\infty}h_R(t')|A(t-t')|^2dt'\right]
```

where:
- ``A(z,t)`` is the complex envelope of the pulse
- ``\alpha(\omega)`` is the linear loss/gain coefficient
- ``\beta_k`` are the dispersion coefficients
- ``\gamma(\omega)`` is the nonlinear coefficient
- ``\omega_0`` is the carrier frequency
- ``f_R`` is the fractional contribution of the Raman effect
- ``h_R(t)`` is the Raman response function

## Physical Effects

### Linear Effects

#### Loss and Gain
The term ``-\frac{\alpha}{2}A`` describes:
- Material absorption
- Scattering losses
- Gain in active media[^Agrawal2019ch6]

#### Dispersion
The dispersion operator ``i\sum_{k\geq 2}\frac{\beta_k}{k!}\left(i\frac{\partial}{\partial t}\right)^k A`` accounts for all orders of chromatic dispersion through the coefficients ``\beta_k``, which come from the Taylor expansion:

```math
\beta(\omega) = \beta_0 + \beta_1(\omega-\omega_0) + \frac{\beta_2}{2}(\omega-\omega_0)^2 + \frac{\beta_3}{6}(\omega-\omega_0)^3 + ...
```

Common values at 1550nm[^Agrawal2019][^Dudley2006]:
- ``\beta_2 \approx -20 \text{ ps}^2/\text{km}`` (anomalous dispersion)
- ``\beta_3 \approx 0.1 \text{ ps}^3/\text{km}``

### Nonlinear Effects

#### Kerr Effect
The instantaneous nonlinear response is characterized by the coefficient ``\gamma``:

```math
\gamma = \frac{2\pi n_2}{\lambda A_\text{eff}}
```

where ``n_2`` is the nonlinear refractive index and ``A_\text{eff}`` is the effective mode area.

Typical values:
- Standard fiber: ``\gamma \approx 1 \text{ W}^{-1}\text{km}^{-1}``[^Agrawal2019]
- Highly nonlinear fiber: ``\gamma \approx 10 \text{ W}^{-1}\text{km}^{-1}``[^Dudley2006]
- Silicon waveguides: ``\gamma > 100 \text{ W}^{-1}\text{m}^{-1}``[^Foster2008]

#### Self-Steepening
The operator ``1 + \frac{i}{\omega_0}\frac{\partial}{\partial t}`` accounts for the intensity dependence of the group velocity[^Agrawal2019ch2]. This becomes important for:
- Ultrashort pulses (< 100 fs)
- High peak powers
- Long propagation distances

#### Raman Effect
The delayed nonlinear response is modeled through:
- Fractional strength ``f_R`` (typically 0.18 in silica)
- Response function ``h_R(t)``

The default implementation uses the Lin-Agrawal model[^LinAgrawal2006]:

```math
h_R(t) = \frac{\tau_1^2 + \tau_2^2}{\tau_1\tau_2^2} \exp(-t/\tau_2)\sin(t/\tau_1)
```

with ``\tau_1 = 12.2 \text{ fs}`` and ``\tau_2 = 32 \text{ fs}`` for silica.

## Mathematical Derivation

For a complete derivation of the GNLSE starting from Maxwell's equations, see the appendix.

## Appendix: GNLSE Derivation

### Starting from Maxwell's Equations

The wave equation in a dielectric medium is:
```math
\nabla^2 \mathbf{E} - \mu_0 \epsilon_0 \frac{\partial^2 \mathbf{E}}{\partial t^2} = \mu_0 \frac{\partial^2 \mathbf{P}}{\partial t^2}
```
where ``\mathbf{P} = \mathbf{P}_{\text{L}} + \mathbf{P}_{\text{NL}}`` includes both linear and nonlinear polarization.

### The Slowly Varying Envelope Approximation

Writing the electric field as:
```math
\mathbf{E}(z, t) = \frac{1}{2} \left[ A(z, t) e^{i(\beta_0 z - \omega_0 t)} + \text{c.c.} \right]
```

and applying the SVEA leads to:
```math
\frac{\partial A}{\partial z} + i \left[ \beta(\omega) - \beta_0 \right] A = i \frac{\omega_0^2}{2 \beta_0 c^2} P_{\text{NL}}
```

### Including Higher-Order Effects

Adding dispersion, Raman effect, and self-steepening leads to the full GNLSE shown above.

[^Agrawal2019]: G. P. Agrawal, "Nonlinear Fiber Optics," Academic Press, 6th Edition (2019).
[^Agrawal2019ch2]: Chapter 2 of Agrawal's "Nonlinear Fiber Optics" covers nonlinear effects in detail.
[^Agrawal2019ch3]: Chapter 3 of Agrawal's "Nonlinear Fiber Optics" covers dispersion in detail.
[^Agrawal2019ch6]: Chapter 6 of Agrawal's "Nonlinear Fiber Optics" covers optical amplifiers.
[^Dudley2006]: J. M. Dudley et al., "Supercontinuum generation in photonic crystal fiber," Rev. Mod. Phys. 78, 1135 (2006).
[^Foster2008]: M. A. Foster et al., "Silicon-chip-based ultrafast optical oscilloscope," Nature 456, 81 (2008).
[^LinAgrawal2006]: Q. Lin and G. P. Agrawal, "Raman response function for silica fibers," Opt. Lett. 31, 3086 (2006).