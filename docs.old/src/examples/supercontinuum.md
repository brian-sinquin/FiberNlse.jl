# Supercontinuum Generation

This example demonstrates supercontinuum generation in a highly nonlinear fiber. A supercontinuum is created when an ultrashort pulse undergoes extreme spectral broadening through various nonlinear processes[^Dudley2006].

```@setup supercontinuum
using FiberNlse
using Plots
using FFTW
include("plotting_setup.jl")
```

## Physical Mechanism

Supercontinuum generation involves multiple nonlinear effects[^Dudley2010]:
1. **Initial stage:**
   - Self-phase modulation
   - Four-wave mixing
2. **Soliton dynamics:**
   - Soliton fission
   - Dispersive wave generation
3. **Late stage:**
   - Raman self-frequency shift
   - Soliton-dispersive wave interaction

## Simulation Setup

```@example supercontinuum
# Grid parameters (from Dudley et al. 2006)
N = 2^13
T = 12.5e-12          # Time window
t = (-N÷2:N÷2-1) * T/N  # Time grid

# Pulse parameters
t₀ = 28.4e-15        # Pulse duration (28.4 fs)
P₀ = 10e3            # Peak power (10 kW)
λ₀ = 835e-9          # Center wavelength

# Create input pulse (ensure ComplexF64)
A = ComplexF64.(@. sqrt(P₀) * sech(t/t₀))

# Fiber parameters (PCF parameters from Dudley et al. 2006)
L = 0.15             # Length (15 cm)
γ = 0.11             # Nonlinear parameter (W⁻¹m⁻¹)
β₂ = -1.1e-26        # GVD (s²/m)
β₃ = 6.9e-41         # Third-order dispersion (s³/m)
β₄ = -9.5e-56        # Fourth-order dispersion (s⁴/m)

# Create waveguide with Raman effect and self-steepening
wg = Waveguide(0.0, [0.0, β₂, β₃, β₄], γ, λ₀, L,
               raman_linagrawaal(),
               true)

# Solve GNLSE
sol = gnlse(A, t, wg, nsaves=500)

# First visualization block: Evolution plots (heatmaps)
p1 = plot_temporal_evolution(sol, title="(a) Temporal Evolution")

p2 = plot_spectral_evolution(sol, title="(b) Spectral Evolution")

# Combine evolution plots
p_evolution = plot(p1, p2,
                  layout=(2,1),
                  size=(800,600),
                  plot_title="Supercontinuum Evolution")

# Second visualization block: Profile plots
z_indices = [1, 100, 250, 500]  # input, early, middle, and final stages
labels = ["Input", "z ≈ 3 cm", "z ≈ 7.5 cm", "Output"]

p3 = plot_temporal_slices(sol, z_indices, labels, title="(a) Temporal Profiles")

p4 = plot_spectral_slices(sol, z_indices, labels, title="(b) Spectral Profiles")

# Combine profile plots
p_profiles = plot(p3, p4,
                 layout=(1,2),
                 size=(900,300),
                 plot_title="Pulse Profiles at Different Distances")

# Display both visualization blocks
p_evolution  # Display evolution plots
p_profiles   # Display profile plots
```

## Analysis

The simulation shows key features of supercontinuum generation[^Dudley2006]:

1. **Temporal domain:**
   - Initial pulse compression (z ≈ 2 cm)
   - Soliton fission (z ≈ 5 cm)
   - Multiple ejected solitons

2. **Spectral domain:**
   - Initial SPM-based broadening
   - Asymmetric broadening due to Raman effect
   - Dispersive wave generation at shorter wavelengths
   - Final bandwidth spanning over 500 nm

The parameters used here match those from the seminal paper by Dudley et al.[^Dudley2006], which has become a standard benchmark for supercontinuum generation simulations.

[^Dudley2006]: J. M. Dudley, G. Genty, and S. Coen, "Supercontinuum generation in photonic crystal fiber," Rev. Mod. Phys. 78, 1135 (2006).
[^Dudley2010]: J. M. Dudley and J. R. Taylor, "Supercontinuum Generation in Optical Fibers," Cambridge University Press (2010).
