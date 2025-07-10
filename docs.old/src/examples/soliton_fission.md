# Soliton Fission

This example demonstrates soliton fission in a highly nonlinear waveguide. Soliton fission is a fundamental process in ultrafast nonlinear optics where higher-order solitons break into their constituent fundamental solitons[^Dudley2006].

```@setup soliton_fission
using FiberNlse
using Plots
include("plotting_setup.jl")
```

## Physical Mechanism

Soliton fission occurs through the following process[^Agrawal2019ch5]:
1. A higher-order soliton (N > 1) forms when input pulse power ``P_0 > P_\text{sol}``, where:
   ```math
   P_\text{sol} = \frac{|\beta_2|}{\gamma T_0^2}
   ```
2. Perturbations (higher-order dispersion, Raman effect) break the periodic evolution
3. The pulse splits into N fundamental solitons with peak powers[^Dudley2006]:
   ```math
   P_k \approx P_0(2N - 2k + 1)^2
   ```
   where k = 1 to N, and the highest power soliton emerges first.

## Simulation Parameters

```@example soliton_fission
# Grid parameters (adapted from Dudley et al. 2006)
N = 2^13
T = 20e-12
t = (-N÷2:N÷2-1) * T/N

# Input pulse parameters
t₀ = 50e-15    # Very short pulse (50 fs)
P₀ = 10e3      # High peak power (10 kW)
A = ComplexF64.(@. sqrt(P₀) * sech(t/t₀))

# Waveguide parameters
L = 0.1        # Length (10 cm)
β₂ = -5e-26    # GVD (anomalous dispersion)
β₃ = 1e-40     # Third order dispersion
γ = 0.1        # Nonlinear coefficient

# Create waveguide with Raman effect and self-steepening
wg = Waveguide(0.0, [0.0, β₂, β₃], γ, 835e-9, L,
               raman_linagrawaal(),
               true)

# Calculate soliton order
N = sqrt(γ * P₀ * t₀^2 / abs(β₂))
println("Soliton order N ≈ $(round(N, digits=1))")

# Solve GNLSE
sol = gnlse(A, t, wg, nsaves=500)

# Create multi-panel visualization
p1 = plot_temporal_evolution(sol, title="(a) Temporal Evolution")

p2 = plot_spectral_evolution(sol, title="(b) Spectral Evolution")

# Plot temporal and spectral slices
z_indices = [1, 100, 250, 500]
labels = ["Input", "z ≈ 3 cm", "z ≈ 7.5 cm", "Output"]
p3 = plot_temporal_slices(sol, z_indices, labels, title="(c) Temporal Profiles")
p4 = plot_spectral_slices(sol, z_indices, labels, title="(d) Spectral Profiles")

# Combine all plots
plot(p1, p2, p3, p4,
     layout=(2,2),
     size=(1000,800),
     plot_title="Soliton Fission (N ≈ $(round(N, digits=1)))")
```

## Analysis

The simulation reveals key features of soliton fission[^Dudley2006][^Gaeta2002]:

1. **Initial Evolution:**
   - Initial temporal compression
   - Spectral broadening due to SPM
   - Maximum compression at ``z_c \approx \pi L_D/2N``, where ``L_D = T_0^2/|\beta_2|``

2. **Fission Process:**
   - Asymmetric breakup due to Raman effect
   - First ejected soliton has highest peak power
   - Each subsequent soliton is weaker and slower

3. **Spectral Features:**
   - Initial symmetric broadening
   - Development of red-shifted Raman solitons
   - Generation of blue-shifted dispersive waves

This process is fundamental to supercontinuum generation in the anomalous dispersion regime[^Dudley2006].

[^Agrawal2019ch5]: Chapter 5 of G. P. Agrawal, "Nonlinear Fiber Optics," 6th ed., Academic Press (2019).
[^Dudley2006]: J. M. Dudley, G. Genty, and S. Coen, "Supercontinuum generation in photonic crystal fiber," Rev. Mod. Phys. 78, 1135 (2006).
[^Gaeta2002]: A. L. Gaeta, "Nonlinear propagation and continuum generation in microstructured optical fibers," Opt. Lett. 27, 924 (2002).
