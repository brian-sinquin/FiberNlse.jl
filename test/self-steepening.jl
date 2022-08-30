include("../src/FiberNlse.jl")
using .FiberNlse
using Plots

# Simulation dimension
Nₜ, Nₗ = (2000, 5000);

# Fiber properties
L = 20e3; # Fiber length

# Signal properties
T = 100e-12; # Signal duration
λ = 1550e-9; # Wavelength
τ = 1e-12; # Pulse duration
N = 1 # Soliton number

fib = FiberNlse.Fiber(
    L,
    FiberNlse.Dispersion([-2.5508963989899877e-26, -0.5508963989899877e-39]),
    1 / 1000,
    0,
    λ,
);
fib.D
t = (0:(Nₜ - 1)) * T / Nₜ .- 0.5T

# Input construction
P₀ = 30#abs((fib.D.β[1]/fib.γ/τ^2)*N^2) # Soliton power
Ψₒ = @. sqrt(P₀) / cosh(t / τ) # Soliton formula

field = FiberNlse.propagate(Ψₒ, [fib, fib], T, Nₗ; progress=true); # run the simulation

begin
    plot(abs2.(field.ψ[end, :]))
    plot!(abs2.(Ψₒ))
end
sum(abs2.(Ψₒ))
sum(abs2.(field.ψ[end, :]))
# Testing soliton propagation (including losses)
@test isapprox(abs2.(Ψₒ .* exp(-fib.α * L)), abs2.(field.ψ[end, :]), atol=1e-5)
