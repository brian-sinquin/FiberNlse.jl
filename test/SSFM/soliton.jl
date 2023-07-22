@testset "Soliton invariance" begin

    # Simulation dimension
    Nₜ, Nₗ = (5000, 500)

    # Fiber properties
    L = 5.0e3 # Fiber length

    # Signal properties
    τ = 40e-12 # Pulse duration
    T = 400e-12 # Signal duration
    λ = 1550e-9 # Wavelength
    N = 1 # Soliton number

    fib = smf28(L, λ)
    fib.α = 0
    t = collect((-Nₜ÷2:Nₜ÷2-1) * T / Nₜ)

    # Input construction
    P₀ =abs((fib.D.β[1] / fib.γ / τ^2) * N^2) # Soliton power
    Ψₒ = sqrt(P₀) * sech.(t ./ τ) # Soliton formula

    field = propagate(Ψₒ, [fib], T, Nₗ) # run the simulation

    # Testing soliton propagation (including losses)
    @test isapprox(abs2.(Ψₒ .* exp(-0.5*fib.α * L)), abs2.(field.ψ[end, :]), atol = 1e-4)
end
