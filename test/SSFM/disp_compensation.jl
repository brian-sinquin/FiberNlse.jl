@testset "Dispersion compensation" begin
    # Simulation dimension
    Nₜ, Nₗ = (2000, 200)

    # Fiber properties
    L = 2.0e3 # Fiber length

    # Signal properties
    T = 100e-12 # Signal duration
    λ = 1550e-9 # Wavelength
    τ = 3e-12 # Pulse duration

    fib1 = Fiber(L, dispersion(18e-6, λ), 0, 0, λ) # Anomalous dispersion
    fib2 = Fiber(L, dispersion(-18e-6, λ), 0, 0, λ) # Normal dispersion

    t = (0:(Nₜ-1)) * T / Nₜ .- 0.5T

    # Input construction
    P₀ = 1e-3
    Ψₒ = @. sqrt(P₀) / cosh(t / τ) # Soliton formula

    field = propagate(Ψₒ, [fib1, fib2], T, Nₗ) # run the simulation

    # Testing signal propagation (including losses)
    @test isapprox(Ψₒ, field.ψ[end, :])
end
