# FiberNlse.jl Documentation

FiberNlse.jl is a Julia package for simulating nonlinear pulse propagation in optical waveguides. It solves the Generalized Nonlinear Schrödinger Equation (GNLSE) using a 4th order Runge-Kutta in the Interaction Picture (ERK4IP) method.

## Applications

The package can simulate pulse propagation in various waveguiding structures:
- Optical fibers (both standard and photonic crystal fibers)
- Integrated optical waveguides
- Other guided-wave optical structures

## Features

- Solve the GNLSE for pulse propagation in optical waveguides
- Support for:
  - Arbitrary-order dispersion
  - Kerr nonlinearity (with optional wavelength dependence)
  - Linear loss or gain (with optional wavelength dependence)
  - Raman effect (with customizable response function)
  - Self-steepening
- Efficient implementation using FFT-based methods
- Built-in visualization tools

## Example Applications

The package includes several examples demonstrating different nonlinear optical phenomena and waveguide types:

- [Soliton Propagation](examples/soliton.md): Fundamental and higher-order soliton dynamics
- [Soliton Fission](examples/soliton_fission.md): Higher-order soliton breakup
- [Pulse Compression](examples/compression.md): Chirped pulse compression
- [Supercontinuum Generation](examples/supercontinuum.md): Extreme spectral broadening
- [Waveguide Examples](examples/waveguides.md): Implementation for different waveguide types

## Getting Started

1. First, [install](installation.md) the package
2. Follow the [Quick Start Guide](quickstart.md) to run your first simulation
3. Read the [Usage Guide](usage.md) for more detailed information
4. Check the [API Reference](api.md) for complete documentation

## Physical Background

For an introduction to the underlying physics and equations, see the [Theory](theory.md) section.
