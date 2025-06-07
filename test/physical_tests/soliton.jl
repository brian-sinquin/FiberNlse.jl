@testset "Soliton Invariance" begin

    # Simulation dimension
    N = 2^13

    # Fiber properties
    L = 5.0e3 # Fiber length

    # Signal properties
    τ = 40e-12 # Pulse duration
    T = 4000e-12 # Signal duration
    λ = 1550e-9 # Wavelength
    n = 1 # Soliton number

    α = 0.0

    # Use keyword arguments for Waveguide constructor
    fib = Waveguide(α=α, β2=-2.6e-26, γ=1.1e-3, λ=λ, L=L)

    t = (-N÷2:N÷2-1) * T / N


    # Input construction
    P₀ = abs((fib.β2 / fib.γ / τ^2) * n^2) # Soliton power
    Ψₒ = sqrt(P₀) * sech.(t ./ τ) .+ 0.0im # Soliton formula

    # Create a Solver object
    solver = Solver(L, 200)

    # Propagate the pulse
    sol = gnlse(Ψₒ, fib, t, solver)

    # Testing soliton propagation
    @test isapprox(abs2.(Ψₒ), abs2.(sol.A[end, :]), atol=1e-4)
end
