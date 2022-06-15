include("../src/FiberNlse.jl")
using .FiberNlse
# Simulation dimension
Nₜ, Nₗ = (1000,1000);

# Fiber properties
L = 5km; # Fiber length
fib2 = FiberNlse.smf28(L)
fib1 = FiberNlse.Fiber(-fib1.D, fib1.α, fib1.γ,L)

# Signal properties
T = 200*ps; # Signal duration
λ = 1550*nm; # Wavelength
τ = 5*ps; # Pulse duration
N = 1 # Soliton number
sim,t,l1 = FiberNlse.configure(Nₜ,Nₗ,fib1, T, λ);
sim2,t,l2 = FiberNlse.configure(Nₜ,10,fib2, T, λ);
# Input construction
Pp = abs((sim.β2/fib1.γ/τ^2)*N^2) # Soliton power
Ψₒ = @. sqrt(Pp)/cosh(t/τ) # Soliton formula
FiberNlse.inputSignal(sim,Ψₒ);
FiberNlse.simulate(sim, true); # run the simulation
FiberNlse.transition(sim,sim2)
FiberNlse.simulate(sim2, true);
# Visualization
l = vcat(l1,l2.+l1[end])
S = vcat(sim.Ψ,sim2.Ψ)
using Plots
gr()
heatmap(l/km, t/ps,abs2.(S'))
