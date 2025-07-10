# Usage Guide

This guide provides an overview of how to use FiberNlse.jl to simulate pulse propagation in optical fibers.

## Core Types

### Waveguide

The `Waveguide` struct describes the propagation conditions in your fiber:

```julia
Waveguide(
    α,      # Linear loss factor (scalar or frequency-dependent array)
    βs,     # GVD orders vector (starting at 2nd order)
    γ,      # Nonlinear factor (scalar or frequency-dependent array)
    λc,     # Center wavelength (m)
    L,      # Physical length (m)
    raman_model = NoRaman,           # Optional: Raman model
    self_steepening = false          # Optional: Enable self-steepening
)
```

### RamanModel

The `RamanModel` struct describes the Raman response of the waveguide:

```julia
RamanModel(
    fr,             # Fractional contribution of Raman effect
    time_response   # Function returning Raman impulse response vs time
)
```

For simulations without Raman effect, use the `NoRaman` constant.

## Solving the GNLSE

The package uses a 4th order Runge-Kutta in the Interaction Picture (ERK4IP) method to solve the Generalized Nonlinear Schrödinger Equation. 

### Basic Usage

Here's a basic example of pulse propagation:

```julia
using FiberNlse

# Define time grid
N = 2^13                   # number of grid points
T = 20e-12                # time window (s)
t = (-N÷2:N÷2-1) * T/N    # time array

# Create input pulse (ensure it's ComplexF64)
A = ComplexF64.(sqrt(1000) * sech.(t ./ 1e-12))

# Define waveguide
wg = Waveguide(0.0, [-20e-27], 0.1, 1550e-9, 1.0)

# Solve GNLSE
sol = gnlse(A, t, wg)
```

### Adding Raman Effect

To include Raman effect:

```julia
# Create Raman model
raman = raman_linagrawaal()  # Using built-in Raman response

# Create waveguide with Raman effect
wg = Waveguide(0.0, [-20e-27], 0.1, 1550e-9, 1.0, raman)
```

### Self-Steepening

To enable self-steepening:

```julia
wg = Waveguide(0.0, [-20e-27], 0.1, 1550e-9, 1.0, NoRaman, true)
```

## Working with Solutions

The `gnlse` function returns a `Solution` struct containing:

- `At`: Matrix of field amplitudes in time domain
- `Af`: Matrix of field amplitudes in frequency domain
- `t`: Time grid
- `z`: Distance grid
- `ω`: Angular frequency grid

Access solution data:

```julia
# Time domain field at input
At_input = sol.At[1,:]

# Time domain field at output
At_output = sol.At[end,:]

# Frequency domain field at any z-position
Af_at_z = sol.Af[z_index,:]
```
