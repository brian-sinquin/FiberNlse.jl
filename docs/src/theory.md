# Derivation of the Nonlinear Schrödinger Equation (NLSE)

The propagation of an optical pulse in a nonlinear and dispersive medium can be described starting from Maxwell's equations. Here, we outline the derivation of the Nonlinear Schrödinger Equation (NLSE) using the Slowly Varying Envelope Approximation (SVEA), following the approach in G. P. Agrawal's *Nonlinear Fiber Optics*.

## 1. Starting from Maxwell's Equations

The wave equation for the electric field ``\mathbf{E}`` in a dielectric medium is:
```math
\nabla^2 \mathbf{E} - \mu_0 \epsilon_0 \frac{\partial^2 \mathbf{E}}{\partial t^2} = \mu_0 \frac{\partial^2 \mathbf{P}}{\partial t^2}
```
where ``\mathbf{P} = \mathbf{P}_{\text{L}} + \mathbf{P}_{\text{NL}}`` is the total polarization, including both linear and nonlinear contributions.

## 2. Linear and Nonlinear Polarization

- **Linear polarization:** ``\mathbf{P}_{\text{L}} = \epsilon_0 \chi^{(1)} \mathbf{E}``
- **Nonlinear polarization:** For Kerr media, ``\mathbf{P}_{\text{NL}} = \epsilon_0 \chi^{(3)} |\mathbf{E}|^2 \mathbf{E}``

## 3. The Slowly Varying Envelope Approximation (SVEA)

Assume the electric field can be written as:
```math
\mathbf{E}(z, t) = \frac{1}{2} \left[ A(z, t) e^{i(\beta_0 z - \omega_0 t)} + \text{c.c.} \right]
```
where ``A(z, t)`` is the slowly varying envelope, ``\omega_0`` is the carrier frequency, and ``\beta_0 = \beta(\omega_0)``.

The polarization can be similarly written as:
```math
\mathbf{P}(z, t) = \frac{1}{2} \left[ P(z, t) e^{i(\beta_0 z - \omega_0 t)} + \text{c.c.} \right]
```
where ``P(z, t)`` is the envelope of the polarization.

## 4. Substituting into the Wave Equation

Substitute the expressions for ``\mathbf{E}(z, t)`` and ``\mathbf{P}(z, t)`` into the wave equation. After some algebra and using the SVEA (neglecting ``\partial^2 A/\partial z^2``), the equation for the envelope becomes:
```math
\left( \frac{\partial}{\partial z} + i \beta_0 \right)^2 A - \frac{1}{c^2} \left( \frac{\partial}{\partial t} - i \omega_0 \right)^2 A = -\mu_0 \frac{\partial^2 P}{\partial t^2} e^{-i(\beta_0 z - \omega_0 t)}
```
Expanding and keeping only first-order derivatives in ``z`` (SVEA), and transforming to the frequency domain, we get:
```math
\frac{\partial A}{\partial z} + i \left[ \beta(\omega) - \beta_0 \right] A = i \frac{\omega_0^2}{2 \beta_0 c^2} P_{\text{NL}}
```
where ``P_{\text{NL}}`` is the nonlinear polarization envelope.

## 5. Including Dispersion

Expand the propagation constant ``\beta(\omega)`` in a Taylor series around ``\omega_0`` to arbitrary order:
```math
\beta(\omega) = \beta_0 + \sum_{n=1}^{\infty} \frac{1}{n!} \beta_n (\omega - \omega_0)^n
```
where
```math
\beta_n = \left. \frac{d^n \beta}{d\omega^n} \right|_{\omega_0}
```
is the ``n``-th order dispersion coefficient.

Transforming back to the time domain, the dispersion operator acts as:
```math
\sum_{n=1}^{\infty} \frac{i^{n-1}}{n!} \beta_n \frac{\partial^n}{\partial t^n} A
```

## 6. Nonlinear Polarization and the Nonlinear Parameter

For a Kerr medium, the nonlinear polarization envelope is:
```math
P_{\text{NL}}(z, t) = \epsilon_0 \chi^{(3)} |A(z, t)|^2 A(z, t)
```
The nonlinear parameter is defined as:
```math
\gamma = \frac{\omega_0 n_2}{c A_{\text{eff}}}
```
where ``n_2`` is the nonlinear refractive index and ``A_{\text{eff}}`` is the effective mode area.

## 7. The Generalized Nonlinear Schrödinger Equation

Combining the above, the envelope ``A(z, t)`` satisfies:
```math
\frac{\partial A}{\partial z} + \sum_{n=1}^{\infty} \frac{i^{n-1}}{n!} \beta_n \frac{\partial^n A}{\partial t^n} = i \gamma |A|^2 A
```
where:
- ``\beta_n`` are the dispersion coefficients,
- ``\gamma`` is the nonlinear parameter.

## 8. Moving to the Retarded Time Frame

In a reference frame moving at the group velocity (``\tau = t - z/v_g``, with ``v_g = 1/\beta_1``), the equation becomes:
```math
\frac{\partial A}{\partial z} + \sum_{n=2}^{\infty} \frac{i^{n-1}}{n!} \beta_n \frac{\partial^n A}{\partial \tau^n} = i \gamma |A|^2 A
```
where the sum starts from ``n=2`` since the first-order term is eliminated in the moving frame.

This is the generalized form of the Nonlinear Schrödinger Equation (NLSE) for pulse propagation in a dispersive and nonlinear medium, including arbitrary high-order dispersion.

## 9. Including the Raman Effect

In real optical fibers, the nonlinear response is not purely instantaneous. The delayed Raman response can be included by modifying the nonlinear term to account for both instantaneous (electronic) and delayed (Raman) contributions. The nonlinear polarization envelope becomes:
```math
P_{\text{NL}}(z, t) = \epsilon_0 \chi^{(3)} \left[ (1 - f_R) |A(z, t)|^2 A(z, t) + f_R A(z, t) \int_{-\infty}^{\infty} h_R(t') |A(z, t - t')|^2 dt' \right]
```
where:
- ``f_R`` is the fractional contribution of the delayed Raman response (typically ``f_R \approx 0.18`` for silica),
- ``h_R(t)`` is the Raman response function.

The generalized NLSE including the Raman effect is:
```math
\frac{\partial A}{\partial z} + \sum_{n=2}^{\infty} \frac{i^{n-1}}{n!} \beta_n \frac{\partial^n A}{\partial \tau^n} = i \gamma \left[ (1 - f_R) |A|^2 A + f_R A \int_{-\infty}^{\infty} h_R(\tau') |A(\tau - \tau')|^2 d\tau' \right]
```
This equation describes pulse propagation in a nonlinear dispersive medium, including arbitrary high-order dispersion and the delayed Raman response.

## 10. Including Self-Steepening

Self-steepening is a higher-order nonlinear effect that arises due to the intensity dependence of the group velocity, leading to an asymmetric spectral broadening of ultrashort pulses. It can be incorporated by adding a term that accounts for the time derivative of the nonlinear response. The generalized NLSE including both Raman effect and self-steepening becomes:
```math
\frac{\partial A}{\partial z} + \sum_{n=2}^{\infty} \frac{i^{n-1}}{n!} \beta_n \frac{\partial^n A}{\partial \tau^n} =
i \gamma \left[ 1 + \frac{i}{\omega_0} \frac{\partial}{\partial \tau} \right] \left( (1 - f_R) |A|^2 A + f_R A \int_{-\infty}^{\infty} h_R(\tau') |A(\tau - \tau')|^2 d\tau' \right)
```
where the operator ``1 + \frac{i}{\omega_0} \frac{\partial}{\partial \tau}`` accounts for self-steepening, and ``\omega_0`` is the carrier frequency.

This equation describes pulse propagation in a nonlinear dispersive medium, including arbitrary high-order dispersion, the delayed Raman response, and self-steepening.

---

**References:**
- G. P. Agrawal, *Nonlinear Fiber Optics*, 5th Edition, Academic Press, 2013.