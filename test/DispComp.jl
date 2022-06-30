include("../src/FiberNlse.jl")
using .FiberNlse
using ForwardDiff
# Simulation dimension
Nₜ, Nₗ = (5000,1000);

# Fiber properties
L = 5km; # Fiber length

fib2 = FiberNlse.smf28(2L)
fib1 = FiberNlse.Fiber(-fib2.D, fib2.α, fib2.γ,L) # perfectly compensating fiber

bundle=FiberNlse.Bundle([fib1,fib2], 1000)

# Signal properties
T = 100*ps; # Signal duration
λ = 1550*nm; # Wavelength
τ = 10*ps; # Pulse duration
N = 1 # Soliton number
sims,t,l = FiberNlse.configure(bundle,Nₜ, T, λ);

# Input construction
Pp = abs((sims[1].β2/fib1.γ/τ^2)*N^2) # Soliton power
Ψₒ = @. sqrt(Pp)/cosh(t/τ) # Soliton formula

S=FiberNlse.simulate(sims,Ψₒ,false);

# Visualization
using Plots
using FFTW
gr()

heatmap(l/km, t/ps,abs2.(S'))

