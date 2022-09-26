include("../src/FiberNlse.jl")
using .FiberNlse


# Simulation dimension
Nₜ, Nₗ = (2000, 2000);

# Fiber properties

D = 0.9ps / nm / km;
α = 0.026 / km;
γ = 10.1 / W / km;
L = 100.0m;

# Signal propoerties
T = 1000ps; # Signal duration
λ = 1550nm; # Wavelength
τ = 5ps; # Pulse duration
N = 1 # Soliton number
sim, t, l = configure(Nₜ, Nₗ, D, γ, α, L, T, λ);
Pp = 0.550;#abs((sim.β2/γ/τ^2)*N^2) # Soliton power
f0 = 1e10;
β = 2
φ = 1.8pi
# Input construction
t = T * 0.5 * range(-1, stop = 1, length = Nₜ); # Time vector
noiseP = sqrt(1e-8) * randn(ComplexF64, length(t))
noisePhi = pi * 1e-14 * randn(ComplexF64, length(t))
Ψₒ = @. noiseP +
   sqrt(Pp) * cos(2pi * f0 * t) * exp(1im * (noisePhi + β * cos(2pi * f0 * t) + φ))



inputSignal(sim, Ψₒ);
simulate(sim, false);

D = 17ps / nm / km;
α = 0.046 / km;
γ = 1.1 / W / km;
L = 5.0km;

sim2, t, l = configure(Nₜ, Nₗ, D, γ, α, L, T, λ);
transition(sim, sim2)
simulate(sim2)


l = range(0, stop = L, length = Nₗ);
using Plots
using FFTW
heatmap((abs2.(sim2.Ψ') .^ 2))
xlabel!("Local time [ps]")

ylabel!("Distance [km]")
fs = 1 / (t[2] - t[1])
f = fftshift(fftfreq(length(t), fs))
plot(f / GHz, abs.(fftshift(fft(sim.Ψ[1, :] .^ 2))), yaxis = :log10)
plot!(f / GHz, abs.(fftshift(fft(sim2.Ψ[end, :] .^ 2))), yaxis = :log10)
#xlims!(-500,500)
hline!(
    [maximum(abs.(fft(sim2.Ψ[end, :] .^ 2))) / 1000],
    linestyle = :dash,
    label = "-30 dB",
)
