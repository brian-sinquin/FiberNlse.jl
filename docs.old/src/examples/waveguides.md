# Waveguide Examples

This section demonstrates how to simulate pulse propagation in different types of waveguiding structures.

## Standard Step-Index Fiber

```@setup waveguides
using FiberNlse
using Plots
include("plotting_setup.jl")  # Common plotting functions and styles
```

```@example waveguides
using FiberNlse
using Plots

# Grid parameters
N = 2^12
T = 20e-12
t = (-N÷2:N÷2-1) * T/N

# Input pulse
P₀ = 100   # Peak power (W)
τ = 1e-12  # Pulse duration (s)
A = ComplexF64.(@. sqrt(P₀) * sech(t/τ))

# Step-index fiber parameters (typical SMF-28)
α = 0.2    # Loss coefficient (dB/km) converted to m⁻¹
L = 1000.0   # Length (m)
β₂ = 20e-27  # GVD at 1550nm (s²/m)
γ = 1.1e-3   # Nonlinear coefficient (W⁻¹m⁻¹)

wg = Waveguide(
    α/4.343,  # Convert from dB/km to m⁻¹
    [0.0, β₂],
    γ,
    1550e-9,  # Operating wavelength
    L
)

sol = gnlse(A, t, wg)

# Plot temporal profiles at input and output
plot_temporal_slices(sol, [1, length(sol.z)], ["Input", "Output"], title="Step-Index Fiber: Temporal Profiles")
```

## Silicon Waveguide

```@example waveguides
# Grid parameters (finer for shorter pulses)
N = 2^13
T = 5e-12
t = (-N÷2:N÷2-1) * T/N

# Input pulse (shorter for silicon photonics)
P₀ = 10    # Peak power (W)
τ = 100e-15  # 100 fs pulse
A = ComplexF64.(@. sqrt(P₀) * sech(t/τ))

# Silicon waveguide parameters
α = 2.0    # Loss coefficient (dB/cm) converted to m⁻¹
L = 0.01   # Length (1 cm)
β₂ = -0.9e-24  # GVD at 1550nm (s²/m)
β₃ = 1e-39    # Third order dispersion (s³/m)
γ = 100.0       # Nonlinear coefficient (W⁻¹m⁻¹)

# Create waveguide with self-steepening (important at these scales)
wg = Waveguide(
    α*100,    # Convert from dB/cm to m⁻¹
    [0.0, β₂, β₃],
    γ,
    1550e-9,
    L,
    NoRaman,  # No Raman in silicon at these timescales
    true      # Enable self-steepening
)

sol = gnlse(A, t, wg)

# Plot temporal and spectral evolution
p1 = plot_temporal_slices(sol, [1, length(sol.z)], ["Input", "Output"], title="(a) Silicon Waveguide: Temporal Profiles")
p2 = plot_spectral_slices(sol, [1, length(sol.z)], ["Input", "Output"], title="(b) Silicon Waveguide: Spectral Profiles")

plot(p1, p2, layout=(2,1), size=(600,800))
```

## Chalcogenide Fiber

```@example waveguides
# Grid parameters
N = 2^13
T = 10e-12
t = (-N÷2:N÷2-1) * T/N

# Input pulse
P₀ = 1000   # Peak power (W)
τ = 200e-15  # 200 fs pulse
A = ComplexF64.(@. sqrt(P₀) * sech(t/τ))

# Chalcogenide fiber parameters (As₂Se₃)
α = 1.0     # Loss (dB/m)
L = 0.1     # Length (10 cm)
β₂ = -400e-27  # Large anomalous GVD
γ = 10.0        # High nonlinearity
fr = 0.03     # Raman fraction (smaller than silica)

# Create waveguide with custom Raman response
function chalcogenide_raman(t)
    τ₁ = 23e-15  # Faster response than silica
    τ₂ = 230e-15
    return @. (τ₁^2 + τ₂^2)/(τ₁*τ₂^2) * exp(-t/τ₂) * sin(t/τ₁) * (t > 0)
end

wg = Waveguide(
    α/4.343,
    [0.0, β₂],
    γ,
    1550e-9,
    L,
    RamanModel(fr, chalcogenide_raman),  # Custom Raman model
    true  # Enable self-steepening
)

sol = gnlse(A, t, wg)

# Plot temporal evolution
plot_temporal_evolution(sol, title="Chalcogenide Fiber: Temporal Evolution")
```

See the [Technical Reference](@ref) section for detailed information about units, parameters, and characteristic lengths used in nonlinear optics simulations.
