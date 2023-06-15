@testset "Soliton invariance" begin
    # Simulation dimension
    Nₜ, Nₗ = (2^10, 400)

    # Fiber properties
    L = 5.0e3 # Fiber length

    # Signal properties
    τ = 10e-12 # Pulse duration
    T = τ*100 # Signal duration
    λ = 1550e-9 # Wavelength
    N = 1 # Soliton number

    fib = smf28(L, λ)
    
    fib.α=0.0
    t = (0:(Nₜ-1)) * T / Nₜ .- 0.5T

    # Input construction
    P₀ = abs((fib.D.β[1] / fib.γ / τ^2) * N^2) # Soliton power
    Ψₒ = @. sqrt(P₀) * sech(t / τ) # Soliton formula

    field = propagate4(Ψₒ, fib, T, Nₗ; progress=true) # run the simulation

    # Testing soliton propagation (including losses)
    @test isapprox(abs2.(Ψₒ), abs2.(field.ψ[end, :]), atol = 1e-4)
   
 end
