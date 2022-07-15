include("../src/FiberNlse.jl")
using .FiberNlse
using Test
# Simulation dimension
Nₜ, Nₗ = (2000,1000);

# Fiber properties
L = 5.0*km; # Fiber length
fib = FiberNlse.smf28(L)

# Signal properties
T = 1000*ps; # Signal duration
λ = 1550*nm; # Wavelength
τ = 100*ps; # Pulse duration
N = 1 # Soliton number
sim,t,l = FiberNlse.configure(Nₜ,Nₗ,fib, T, λ);
# Input construction
P₀ =  abs((sim.β2/fib.γ/τ^2)*N^2) # Soliton power
Ψₒ = @. sqrt(P₀)/cosh(t/τ) # Soliton formula

FiberNlse.inputSignal(sim,Ψₒ);

FiberNlse.simulate(sim, false); # run the simulation

# Testing soliton propagation (including losses)
@test isapprox(abs2.(Ψₒ.*exp(-0.5*fib.α*L)), abs2.(sim.Ψ[end,:]), atol=0.005)

@time sol=FiberNlse.simulate2(sim);
sol=sol(l)
# Visualization

#using Plots
#heatmap(t/ps, l/km, abs2.(sim.Ψ))
#surface(l/km,t/ps, abs2.(sim.Ψ'), camera=(60,20))
#plot(abs2.(sim.Ψ[end,:]))
#plot!(abs2.(sol[:,end]))