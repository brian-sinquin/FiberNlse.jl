# Quick Start Guide

This guide will help you get started with FiberNlse.jl quickly. We'll simulate a simple soliton propagation in a nonlinear fiber.

## Basic Example

```@example quickstart
using FiberNlse
using Plots

# 1. Define simulation parameters
N = 2^12                   # number of grid points
T = 20e-12                # time window width [s]
t = (-N÷2:N÷2-1) * T/N   # time grid

# 2. Create input pulse (sech pulse)
power = 1000              # peak power [W]
t0 = 1e-12               # pulse duration [s]
A = ComplexF64.(@. sqrt(power) * sech(t / t0))  # Ensure ComplexF64 type

# 3. Set fiber parameters
length = 0.1              # fiber length [m]
β₂ = -20e-27             # GVD parameter [s²/m]
γ = 0.1                  # nonlinear parameter [1/W/m]

# 4. Create waveguide
wg = Waveguide(0.0,      # linear loss
               [0.0, β₂], # dispersion terms (starting from β₂)
               γ,         # nonlinearity
               1550e-9,   # center wavelength
               length)    # fiber length

# 5. Solve GNLSE
sol = gnlse(A, t, wg)

# 6. Plot results
plot(t*1e12, abs2.(sol.At[1,:]), 
     xlabel="Time (ps)", 
     ylabel="Power (W)",
     label="Input")
plot!(t*1e12, abs2.(sol.At[end,:]),
      label="Output")
```

## Understanding the Results

The above example simulates the propagation of a sech pulse through a fiber with:
- Anomalous dispersion (β₂ < 0)
- Kerr nonlinearity (γ ≠ 0)
- No linear loss (α = 0)
- No Raman effect (default)
- No self-steepening (default)

For these parameters, the pulse forms a fundamental soliton, maintaining its shape during propagation because:
1. The peak power and pulse duration are chosen to balance dispersion and nonlinearity
2. The fundamental soliton condition N = 1 is approximately satisfied

## Adding Advanced Effects

To include Raman effect and self-steepening:

```@example quickstart
# Create waveguide with Raman effect and self-steepening
wg_advanced = Waveguide(0.0, [0.0, β₂], γ, 1550e-9, length,
                       raman_linagrawaal(),  # Add Raman
                       true)  # Enable self-steepening

# Solve and plot
sol_adv = gnlse(A, t, wg_advanced)

plot(t*1e12, abs2.(sol_adv.At[end,:]), 
     xlabel="Time (ps)", 
     ylabel="Power (W)",
     label="With Raman + Self-steepening")
```

## Next Steps

To learn more:
- See the [Usage Guide](usage.md) for detailed API information
- Check the [Theory](theory.md) section for physical background
- Try the example scripts in the `examples/` directory for more complex scenarios
- Learn about [visualization techniques](visualization.md) for analyzing results
