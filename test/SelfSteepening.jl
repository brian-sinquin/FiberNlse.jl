include("../src/FiberNlse.jl")
using Plots;
gr()

U = FiberNlse # Alias for units;

# Simulation dimension
Nₜ, Nₗ = (1000,1000);

# Fiber properties

D = -17*U.ps/U.nm/U.km;
α = 0.046/U.km;
γ =  1.1/U.W/U.km;
L = 5.0*U.km;
N=1
fib = U.Fiber(L, α, D, γ)


# Signal propoerties
T = 100*U.ps; # Signal duration
λ = 1550*U.nm; # Wavelength
τ = 5 * U.ps; # Pulse duration
sim, t = U.configure(Nₜ,Nₗ,D, γ, α, L, T, λ);
Pp = abs((sim.β2/γ/τ^2)*N^2) # Soliton power

# Input construction
Ψₒ = Vector{ComplexF64}(@. sqrt(Pp)/cosh(t/τ)) # Soliton formula





U.inputSignal(sim,Ψₒ);
U.simulate(sim, false);

l = range(0,stop=L, length=Nₗ);

heatmap(abs.(sim.Ψ).^2)
xlabel!("Local time [ps]")

ylabel!("Distance [km]")
