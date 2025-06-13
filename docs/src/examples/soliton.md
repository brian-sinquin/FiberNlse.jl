# Optical Soliton Propagation

This example demonstrates the propagation of optical solitons - remarkable waves that maintain their shape during propagation due to a perfect balance between dispersion and nonlinearity. First discovered by Hasegawa and Tappert in 1973[^1], solitons have revolutionized our understanding of nonlinear wave phenomena and found applications in optical communications and ultrafast lasers.

```@setup soliton
using FiberNlse
using Plots
include("plotting_setup.jl")  # Common plotting functions and styles
```

## Mathematical Foundation

The normalized nonlinear Schrödinger equation (NLSE) governing pulse propagation in optical fibers is:

```math
i\frac{\partial \psi}{\partial \xi} + \frac{1}{2}\frac{\partial^2 \psi}{\partial \tau^2} + |\psi|^2\psi = 0
```

where ``\xi = z/L_D`` is the normalized distance (``L_D = T_0^2/|\beta_2|`` is the dispersion length) and ``\tau = t/T_0`` is the normalized time.

The fundamental soliton solution (N=1) is:

```math
\psi(\xi, \tau) = \text{sech}(\tau)e^{i\xi/2}
```

In physical units, a soliton of order N is described by:

```math
A(z, t) = N \sqrt{P_0} \text{sech}\left(\frac{t}{\tau}\right)\exp\left(i\frac{z}{2L_D}\right)
```

where:
- ``P_0 = \frac{| \beta_2 |}{\gamma \tau^2}`` is the peak power of the soliton,
- ``\tau`` is the pulse duration, which determines the temporal width of the soliton,
- ``\gamma`` is the nonlinear coefficient, which quantifies the strength of the nonlinearity,
- ``\beta_2`` is the group velocity dispersion (GVD) parameter, which quantifies the dispersion in the fiber,
- ``z`` is the propagation distance,
- ``t`` is the time,
- ``N`` is the soliton order, which determines the balance between dispersion and nonlinearity.

The soliton order ``N`` is given by:

```math
N = \sqrt{\frac{\gamma P_0 \tau^2}{|\beta_2|}}
```

For a fundamental soliton, ``N = 1``, meaning the dispersion and nonlinearity are perfectly balanced, resulting in a stable pulse that does not change shape during propagation. Higher-order solitons (``N > 1``) exhibit periodic evolution, where the pulse periodically compresses and broadens as it propagates.

Solitons are not limited to optics. They are found in various fields of science, including fluid dynamics (e.g., water waves), plasma physics, and even biology (e.g., nerve pulse propagation). These phenomena share the common feature of maintaining their shape due to a balance between competing effects.

*References:*
1. P. G. Agrawal, *Nonlinear Fiber Optics*[^1].
2. A. Hasegawa and F. Tappert, "Transmission of stationary nonlinear optical pulses in dispersive dielectric fibers. I. Anomalous dispersion," *Applied Physics Letters*[^2].
3. Y. S. Kivshar and G. P. Agrawal, *Optical Solitons: From Fibers to Photonic Crystals*[^3].

[^1]: P. G. Agrawal, *Nonlinear Fiber Optics*, Academic Press, 5th Edition, 2012.
[^2]: A. Hasegawa and F. Tappert, "Transmission of stationary nonlinear optical pulses in dispersive dielectric fibers. I. Anomalous dispersion," *Applied Physics Letters*, 23(3), 1973.
[^3]: Y. S. Kivshar and G. P. Agrawal, *Optical Solitons: From Fibers to Photonic Crystals*, Academic Press, 2003.

## Simulation Parameters

Define the simulation parameters, including the fiber properties and signal properties. These parameters are used to model the propagation of solitons in the fiber.

```@example soliton
# Simulation dimension
N = 2^13

# Fiber properties
L = 30.0e3 # Fiber length in meters

# Signal properties
τ = 15e-12 # Pulse duration in seconds
T = 500e-12 # Signal duration in seconds
λ = 1550e-9 # Wavelength in meters

α = 0.0 # Fiber loss (assumed to be zero for simplicity)

fib = Waveguide(α, [0.0, -2.6e-26], 1.1e-3, λ, L) # Define the waveguide
t = (-N÷2:N÷2-1) * T / N # Time grid
```

## Solve the GNLSE

Propagate the pulse through the fiber using the generalized nonlinear Schrödinger equation (GNLSE) for soliton orders 1, 2, and 3. The GNLSE models the interplay between dispersion, nonlinearity, and other effects in the fiber (only dispersion and SPM here).

```@example soliton
# Order 1
P₀_1 = abs((fib.βs[2] / fib.γ / τ^2) * 1^2) # Soliton power for order 1
# Ensure the input pulse is a ComplexF64 vector
Ψₒ_1 = ComplexF64.(@. sqrt(P₀_1) * sech(t ./ τ))
sol_1 = gnlse(Ψₒ_1, t, fib, nsaves=200)

# Order 2
P₀_2 = abs((fib.βs[2] / fib.γ / τ^2) * 2^2) # Soliton power for order 2
Ψₒ_2 = sqrt(P₀_2) * sech.(t ./ τ) .+ 0.0im # Soliton formula for order 2
sol_2 = gnlse(Ψₒ_2, t, fib, nsaves=200)

# Order 3
P₀_3 = abs((fib.βs[2] / fib.γ / τ^2) * 3^2) # Soliton power for order 3
Ψₒ_3 = sqrt(P₀_3) * sech.(t ./ τ) .+ 0.0im # Soliton formula for order 3
sol_3 = gnlse(Ψₒ_3, t, fib, nsaves=200)
```

## Plotting the Results

```@example soliton
# First visualization block: Temporal evolution heatmaps
p1 = plot_temporal_evolution(sol_1, title="(a) Order 1 Soliton")
plot!(p1, xlims=(-50, 50))

p2 = plot_temporal_evolution(sol_2, title="(b) Order 2 Soliton")
plot!(p2, xlims=(-50, 50))

p3 = plot_temporal_evolution(sol_3, title="(c) Order 3 Soliton")
plot!(p3, xlims=(-50, 50))

# Combine heatmaps
p_heat = plot(p1, p2, p3,
             layout=(1,3),
             size=(900,300),
             plot_title="Temporal Evolution of Different Order Solitons")

# Second visualization block: Temporal slices
z_indices = [1, length(sol_1.z)÷4, length(sol_1.z)÷2]
labels = ["z = 0", "z = L/4", "z = L/2"]

p4 = plot_temporal_slices(sol_1, z_indices, labels, title="(a) Order 1")
plot!(p4, xlims=(-50, 50))

p5 = plot_temporal_slices(sol_2, z_indices, labels, title="(b) Order 2")
plot!(p5, xlims=(-50, 50))

p6 = plot_temporal_slices(sol_3, z_indices, labels, title="(c) Order 3")
plot!(p6, xlims=(-50, 50))

# Combine line plots
p_slices = plot(p4, p5, p6,
                layout=(1,3),
                size=(900,300),
                plot_title="Temporal Profiles at Different Distances")

# Display both plots
p_heat  # Display heatmaps
p_slices  # Display line plots
```

## Results and Discussion

- **Order 1 (a):** The fundamental soliton (``N = 1``) maintains its shape during propagation. This is expected as the balance between dispersion and nonlinearity is perfect, resulting in a stable pulse that does not change over the fiber length.

- **Order 2 (b):** The second-order soliton (``N = 2``) exhibits periodic evolution. The pulse compresses and broadens as it propagates, demonstrating the characteristic periodic behavior of higher-order solitons.

- **Order 3 (c):** The third-order soliton (``N = 3``) shows even more pronounced periodic evolution compared to the second-order soliton. The pulse undergoes stronger compression and broadening cycles, highlighting the increasing complexity of higher-order soliton dynamics.

These results illustrate the rich dynamics of solitons in optical fibers and the importance of accurately modeling the interplay between dispersion and nonlinearity. The GNLSE provides a powerful framework for understanding these phenomena.

