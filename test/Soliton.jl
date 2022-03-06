include("../src/FiberNlse.jl")
using Plots;
U = NLSE # Alias for units;

# Simulation dimension
Nₜ, Nₗ = (1000,1000);

# Fiber properties

D = -17*U.ps/U.nm/U.km;
α = 0.046/U.km;
γ =  1.1/U.W/U.km;
L = 5.0*U.km;

# Signal propoerties
T = 20*U.ps; # Signal duration
λ = 1550*U.nm; # Wavelength
τ = 2 * U.ps # Pulse duration
N= 3 # Soliton number
Pp = abs((sim.β2/γ/τ^2)*N^2) # Soliton power

# Input construction
t = T*0.5*range(-1, stop=1, length=Nₜ) # Time vector
Ψₒ = Vector{ComplexF32}(@. 0.5*sqrt(Pp)/cosh(t/τ)) # Soliton formula


sim = U.configure(Nₜ,Nₗ,D, γ, α, L, T, λ)


U.initialSignal(sim,Ψₒ)
U.simulate(sim)

l = range(0,stop=L, length=Nₗ)

heatmap(abs.(sim.Ψ).^2)
xlabel!("Local time [ps]")
ylabel!("Distance [km]")


