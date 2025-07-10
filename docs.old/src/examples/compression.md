# Pulse Compression

This example demonstrates temporal pulse compression in a nonlinear waveguide, a fundamental technique in ultrafast optics. The compression is achieved through the interplay between anomalous dispersion and nonlinear phase modulation of a chirped pulse.

```@setup compression
using FiberNlse
using Plots
include("plotting_setup.jl")  # Common plotting functions and styles
```

## Theory

The pulse compression process relies on the combined effects of group velocity dispersion (GVD) and self-phase modulation (SPM). For a pulse with an initial electric field envelope ``A(0,t)``, the evolution is governed by the nonlinear Schrödinger equation:

```math
\frac{\partial A}{\partial z} = -i\frac{\beta_2}{2}\frac{\partial^2 A}{\partial t^2} + i\gamma|A|^2A
```

where ``\beta_2`` is the GVD parameter and ``\gamma`` is the nonlinear coefficient.

The initial chirped Gaussian pulse is described by:

```math
A(0,t) = \sqrt{P_0}\exp\left(-\frac{t^2}{2t_0^2}(1 + iC)\right)
```

where:
- ``P_0`` is the peak power
- ``t_0`` is the initial pulse duration
- ``C`` is the chirp parameter

The compression mechanism works through three stages:
1. Initial positive chirp creates a frequency-time correlation (frequency increases with time)
2. Anomalous dispersion (``\beta_2 < 0``) causes higher frequencies to travel slower
3. SPM induces additional intensity-dependent phase shift, enhancing compression

The optimal compression distance ``z_{opt}`` can be estimated as[^1]:

```math
z_{opt} \approx \frac{t_0^2}{|\beta_2|}\frac{C}{1 + C^2}
```

[^1]: G. P. Agrawal, "Nonlinear Fiber Optics," Academic Press, 6th Edition (2019).
[^2]: W. J. Tomlinson et al., "Compression of optical pulses chirped by self-phase modulation in fibers," JOSA B, Vol. 1, No. 2 (1984).

## Simulation Setup

```@example compression
# Grid parameters
N = 2^13
T = 10e-12
t = (-N÷2:N÷2-1) * T/N

# Input pulse parameters
t₀ = 1e-12     # Initial duration (1 ps)
P₀ = 5000      # Peak power (5 kW)
C = 3.0        # Chirp parameter (positive)

# Create chirped gaussian pulse
A = ComplexF64.(@. sqrt(P₀) * exp(-0.5(t/t₀)^2 * (1 + im*C)))

# Waveguide parameters
L = 0.5        # Length (50 cm)
β₂ = -20e-27   # GVD (anomalous dispersion)
γ = 0.1        # Nonlinear coefficient

# Create waveguide (no Raman or self-steepening needed for compression)
wg = Waveguide(0.0, [0.0, β₂], γ, 1550e-9, L)

# Solve GNLSE with many save points to see evolution
sol = gnlse(A, t, wg, nsaves=500)

# Calculate compression parameters
peak_power = maximum(abs2.(sol.At), dims=2)[:]
comp_factor = peak_power ./ P₀
z_opt = sol.z[argmax(comp_factor)] * 100  # in cm
max_comp = maximum(comp_factor)

println("Maximum compression factor: $(round(max_comp, digits=1))")
println("Optimal distance: $(round(z_opt, digits=1)) cm")

# First visualization block: Evolution plot
p1 = plot_temporal_evolution(sol, title="Temporal Evolution")

# Calculate FWHM at optimal compression
function calculate_fwhm(t, power)
    max_power = maximum(power)
    above_half = power .>= max_power/2
    t_above = t[above_half]
    return (maximum(t_above) - minimum(t_above))
end

fwhm_in = calculate_fwhm(t, abs2.(sol.At[1,:]))
fwhm_opt = calculate_fwhm(t, abs2.(sol.At[argmax(comp_factor),:]))
compression_ratio = fwhm_in/fwhm_opt

println("FWHM Compression ratio: $(round(compression_ratio, digits=1))")

# Display evolution plot
plot!(p1, plot_title="Pulse Evolution During Compression")

# Second visualization block: Analysis plots
# Compression factor plot
p2 = plot(sol.z*100, comp_factor,
          xlabel="Distance (cm)",
          ylabel="Compression Factor",
          label="")
plot!(p2, title="(a) Compression Factor")
vline!([z_opt], label="Optimal Length", linestyle=:dash)

# Temporal profiles at key distances
z_indices = [1,                              # Input
            argmax(comp_factor),             # Optimal compression
            length(sol.z)]                   # Output
labels = ["Input", "Optimal Compression", "Output"]

p3 = plot_temporal_slices(sol, z_indices, labels, title="(b) Temporal Profiles")

# Combine analysis plots
p_analysis = plot(p2, p3,
                 layout=(1,2),
                 size=(900,300),
                 plot_title="Compression Analysis ($(round(compression_ratio, digits=1))× FWHM reduction)")

# Display both blocks
p1         # Display evolution plot
p_analysis # Display analysis plots
```

## Analysis

The simulation demonstrates key features of pulse compression:

1. **Temporal Evolution (a):**
   - Initial pulse shortening
   - Optimal compression at a specific distance
   - Pulse breakup beyond the optimal point

2. **Compression Factor (b):**
   - Rapid increase in peak power
   - Clear optimal compression distance
   - Degradation after optimal point

3. **Temporal Profiles (c):**
   - Input chirped pulse
   - Maximally compressed pulse
   - Post-compression evolution

The compression achieved here is typical for chirped pulse compression in optical fibers, with:
- Significant peak power enhancement
- Substantial pulse duration reduction
- Clean compressed pulse shape at optimal distance
