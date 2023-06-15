@testset "Self-Phase Modulation" begin

    FiberNlse.setPhaseConvention(:positive)
    # Simulation dimension
    Nz = 300

    # Fiber properties
    γ = 1e-3
    L = 100.0

    # Signal properties
    λ = 1550e-9 # Wavelength
    Pp = 1

    T = 100e-12
    Nt = 2^10
    t = T * (0:(Nt-1)) / Nt .- 0.5T# Time vector
    t0 = 10e-12
    ψ₀ = @.  sqrt(Pp) * (1e-4 .+ exp(-(t/t0)^2))

    fib = Fiber(L, dispersion(0.0, λ), γ, 0.0, λ)
    field = propagate(ψ₀, fib, T, Nz)

    φ = phase(field.ψ[end, :])
    φₜₕ =  -γ * abs2.(ψ₀) * L

     @test isapprox(φₜₕ,φ, rtol=1e-4)

end


