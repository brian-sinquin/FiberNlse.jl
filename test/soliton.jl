@testset "Soliton invariance" begin

    # Simulation dimension
    Nₜ, Nₗ = (2000,5000);

    # Fiber properties
    L = 5.0e3; # Fiber length

    # Signal properties
    T = 1000e-12; # Signal duration
    λ = 1550e-9; # Wavelength
    τ = 100e-12; # Pulse duration
    N = 1 # Soliton number

    fib = smf28(L, λ)
    t = (0:Nₜ-1)*T/Nₜ .- 0.5T

    # Input construction
    P₀ =  abs((fib.D.β[1]/fib.γ/τ^2)*N^2) # Soliton power
    Ψₒ = @. sqrt(P₀)/cosh(t/τ) # Soliton formula



    field=propagate(Ψₒ , [fib,fib], T, Nₗ) # run the simulation

    # Testing soliton propagation (including losses)
    @test isapprox(abs2.(Ψₒ.*exp(-fib.α*L)), abs2.(field.ψ[end,:]), atol=1e-5)

end
