# Usage

This section provides a tutorial on how to use the FiberNlse.jl API to set up and run simulations.

## Installation

First, install the FiberNlse.jl package:

```julia
using Pkg
Pkg.add("FiberNlse")
```

## Setting up a simulation

A typical simulation involves the following steps:

1.  **Define the simulation parameters:** This includes the time window, number of points, and propagation step size.
2.  **Create the initial pulse:** Define the initial field profile.
3.  **Define the waveguide:** Specify the fiber parameters, such as dispersion, nonlinearity, and loss.
4.  **Create a solver:** Set up the numerical solver with appropriate parameters.
5.  **Run the simulation:** Propagate the pulse through the fiber using the GNLSE solver.
6.  **Analyze the results:** Extract and visualize the simulation output.

## Example: Simulating Supercontinuum Generation

Here's a step-by-step example demonstrating supercontinuum generation in an optical fiber:

### 1. Define Simulation Parameters

```julia
using FiberNlse

# Simulation parameters
dz = 0.01
nz = 200
nt = 2048
T = 20.0

# Time grid
t = range(-T / 2, T / 2; length=nt)
```

### 2. Create the Initial Pulse

```julia
# Initial pulse (sech-shaped)
A0 = Base.sech.(t)
```

### 3. Define the Waveguide

```julia
# Waveguide parameters
β2 = -0.5
β3 = 0.1
γ = 2.0

# Create the waveguide
wg = WaveGuide(β2=β2, β3=β3, γ=γ)
```

### 4. Create a Solver

```julia
# Create the solver
solver = Solver(dz, nz)
```

### 5. Run the Simulation

```julia
# Propagate the pulse
sol = gnlse(A0, wg, t, solver)
```

### 6. Analyze the Results

```julia
# Access the results
A = sol.A  # Field evolution
z = sol.z  # Propagation distances
```

### 7. Visualize the Results (Optional)

To visualize the results, you can use the `Plots` package:

```julia
using Plots
using FFTW

# Calculate the spectrum at the end of the fiber
ω = fftshift(fftfreq(nt, t[2]-t[1]))
spectrum = abs.(fftshift(fft(A[:, end])))

# Plot the spectrum
plot(ω, spectrum, xlabel="Frequency", ylabel="Intensity", title="Output Spectrum")
```

This example provides a basic framework for setting up and running simulations with FiberNlse.jl. You can modify the parameters, pulse shapes, and waveguide properties to explore different nonlinear effects in optical fibers. Refer to the API documentation for more details on the available functions and options.
