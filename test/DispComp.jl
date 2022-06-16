include("../src/FiberNlse.jl")
using .FiberNlse
# Simulation dimension
Nₜ, Nₗ = (1000,1000);

# Fiber properties
L = 5km; # Fiber length
fib2 = FiberNlse.smf28(L)
fib1 = FiberNlse.Fiber(0, fib2.α, fib2.γ,L)

b=FiberNlse.Bundle([FiberNlse.edfa(15,15)], 1000)
# Signal properties
T = 200*ps; # Signal duration
λ = 1550*nm; # Wavelength
τ = 5*ps; # Pulse duration
N = 1 # Soliton number
sims,t,l = FiberNlse.configure(b,Nₜ, T, λ);
# Input construction
Pp = 1e-3#abs((sims[1].β2/fib1.γ/τ^2)*N^2) # Soliton power
Ψₒ = @. sqrt(Pp)/cosh(t/τ) # Soliton formula
S=FiberNlse.simulate(sims,Ψₒ,false);

# Visualization
using Plots
gr()
heatmap(l/km, t/ps,abs2.(S'))
