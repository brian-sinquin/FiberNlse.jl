using Plots
using FiberNlse
# Simulation dimension
Nₜ, Nₗ = (2000, 300);

# Fiber properties
L = 10e3; # Fiber length

# Signal properties
T = 1000e-12; # Signal duration
λ = 1550e-9; # Wavelength
τ = 15e-12; # Pulse duration
N = 1 # Soliton number

fib = FiberNlse.Fiber(
    L,
    FiberNlse.Dispersion([-10.5508963989899877e-26, 0*-0.5508963989899877e-39]),
    0*1 / 1000,
    0,
    λ,
);
fib.D
t = (0:(Nₜ-1)) * T / Nₜ .- 0.5T

# Input construction
P₀ = 30#abs((fib.D.β[1]/fib.γ/τ^2)*N^2) # Soliton power
Ψₒ = @. sqrt(P₀) / cosh(t / τ) # Soliton formula

field = FiberNlse.propagate(Ψₒ, [fib, fib], T, Nₗ; progress = true); # run the simulation

begin
    plot(real.(field.ψ[end, :]))
    plot!(real.(Ψₒ))
end
sum(real.(Ψₒ))
sum(real.(field.ψ[end, :]))
# Testing soliton propagation (including losses)
@test isapprox(abs2.(Ψₒ .* exp(-fib.α * L)), abs2.(field.ψ[end, :]), atol = 1e-5)
