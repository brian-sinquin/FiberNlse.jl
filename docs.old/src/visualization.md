# Visualization Guide

This guide covers various techniques for visualizing simulation results from FiberNlse.jl.

## Time Domain Plots

### Pulse Evolution
```@example visualization
using Plots
using FiberNlse

# Define simulation parameters
N = 2^13
T = 500e-12
τ = 15e-12
λ = 1550e-9
α = 0.0
fib = Waveguide(α, [0.0, -2.6e-26], 1.1e-3, λ, 30.0e3)
t = (-N÷2:N÷2-1) * T / N

# Generate input and output pulses
Ψₒ = sqrt(abs(fib.βs[2] / fib.γ / τ^2)) * sech.(t ./ τ) .+ 0.0im
sol = gnlse(Ψₒ, t, fib, nsaves=200)

# Plot input and output pulse intensity vs time
plot(t .* 1e12, abs2.(sol.At[1,:]),
     xlabel="Time (ps)",
     ylabel="Power (W)",
     label="Input")
plot!(t .* 1e12, abs2.(sol.At[end,:]),
      label="Output")
```

## Frequency Domain Plots

### Spectrum Evolution
```@example visualization
using FFTW

# Calculate frequency grid
ω = fftshift(fftfreq(length(t), 1 / (t[2] - t[1]))) .* 2π
S_input = abs2.(fftshift(fft(sol.At[1,:])))
S_output = abs2.(fftshift(fft(sol.At[end,:])))

# Plot spectrum in dB scale
plot(ω ./ 2π * 1e-12, 10 * log10.(S_input .+ eps()),
     xlabel="Frequency (THz)",
     ylabel="Spectral Power (dB)",
     label="Input")
plot!(ω ./ 2π * 1e-12, 10 * log10.(S_output .+ eps()),
      label="Output")
```

## Propagation Heatmaps

### Heatmap of Pulse Propagation
```@example visualization
# Generate propagation heatmap
heatmap(t .* 1e12, sol.z, abs2.(sol.At),
        xlabel="Time (ps)",
        ylabel="Distance (m)",
        title="Propagation Heatmap",
        color=:viridis)
```

