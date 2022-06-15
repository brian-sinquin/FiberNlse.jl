include("../src/FiberNlse.jl")
using .FiberNlse


# Simulation dimension
Nₜ, Nₗ = (2000,2000);

# Fiber properties

D = -0.9*U.ps/U.nm/U.km;
α = 0.026/U.km;
γ =  10.1/U.W/U.km;
L = 100.0*U.m;

# Signal propoerties
T = 1000*U.ps; # Signal duration
λ = 1550*U.nm; # Wavelength
τ = 5 * U.ps; # Pulse duration
N = 1 # Soliton number
sim = U.configure(Nₜ,Nₗ,D, γ, α, L, T, λ);
Pp = 0.550;#abs((sim.β2/γ/τ^2)*N^2) # Soliton power
f0=1e10;
β = 2
φ = 1.8pi
# Input construction
t = T*0.5*range(-1, stop=1, length=Nₜ); # Time vector
noiseP = sqrt(1e-8)*randn(ComplexF64,length(t))
noisePhi = pi*1e-14*randn(ComplexF64,length(t))
Ψₒ = noiseP + sqrt(Pp)*cos(2pi*f0*t)*exp(1im*(noisePhi+β*cos(2pi*f0*t)+φ))



U.inputSignal(sim,Ψₒ);
U.simulate(sim, false);

D = -17*U.ps/U.nm/U.km;
α = 0.046/U.km;
γ =  1.1/U.W/U.km;
L = 5.0*U.km;

sim2 = U.configure(Nₜ,Nₗ,D, γ, α, L, T, λ);
U.transition(sim, sim2)



l = range(0,stop=L, length=Nₗ);

heatmap(abs.(sim2.Ψ).^2)
xlabel!("Local time [ps]")

ylabel!("Distance [km]")
fs = 1/(t[2]-t[1])
f = fftshift(fftfreq(length(t), fs))
plot(f/U.GHz,abs.(fftshift(fft(sim.Ψ[1,:].^2))), yaxis=:log10)
plot!(f/U.GHz,abs.(fftshift(fft(sim2.Ψ[end,:].^2))),yaxis=:log10)
#xlims!(-500,500)
hline!([maximum(abs.(fft(sim2.Ψ[end,:].^2)))/1000], linestyle=:dash, label="-30 dB")