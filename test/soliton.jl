@testset "Soliton invariance" begin

    # Simulation dimension
    Nₜ, Nₗ = (5000, 1000)

    # Fiber properties
    L = 5.0e3 # Fiber length

    # Signal properties
    τ = 20e-12 # Pulse duration
    T = Nₜ*τ*10 # Signal duration
    λ = 1550e-9 # Wavelength
    N = 1 # Soliton number

    fib = smf28(L, λ)
    fib.α = 0
    t = (0:(Nₜ-1)) * T / Nₜ .- 0.5T

    # Input construction
    P₀ = abs((fib.D.β[1] / fib.γ / τ^2) * N^2) # Soliton power
    Ψₒ = @. sqrt(P₀) / cosh(t / τ) # Soliton formula

    field = propagate(Ψₒ, [fib], T, Nₗ) # run the simulation

    # Testing soliton propagation (including losses)
    @test isapprox(abs2.(Ψₒ .* exp(-0.5*fib.α * L)), abs2.(field.ψ[end, :]), atol = 1e-5)
end
