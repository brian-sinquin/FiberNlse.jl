using FiberNlse
using Plots
Nₜ, Nₗ = (2^13, 400)

# Fiber properties
L = 50000.0 # Fiber length

# Signal properties
τ = 2e-12 # Pulse duration
T = τ*100 # Signal duration
λ = 1550e-9 # Wavelength
N = 1 # Soliton number

fib = smf28(L, λ)



t = (0:(Nₜ-1)) * T / Nₜ .- 0.5T

# Input construction
P₀ = 1# Soliton power
Ψₒ =  sqrt(P₀) .* (1.0   .+ 1e-5randn(ComplexF64,length(t))) # Soliton formula

field = propagate4(Ψₒ, fib, T, Nₗ; progress=true) # run the simulation

plot(Ψₒ .|> abs2)
plot!(output(field) .|> abs2)

heatmap(abs2.(field.ψ))