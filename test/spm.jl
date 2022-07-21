@testset "Self-Phase Modulation" begin
using DSP

# Simulation dimension
Nₜ, Nₗ = (2049,5000);

# Fiber properties

D = 0*0.9e-6;
α = 0*0.026e-3;
γ =  10.1e-3;
L = 100.0;

# Signal propoerties
λ = 1550e-9; # Wavelength
Pp = 0.550
f₀ = 10.0e9
T = 2/(f₀) # Signal duration




t = T*(0:Nₜ-1)/Nₜ # Time vector
ψ₀ = @. 0*1im .+ sqrt(Pp)*cos(2pi*f₀*t)

fib = Fiber(L, dispersion(D, λ),γ,α,λ)
field = propagate(ψ₀,fib,T,Nₗ)

φ = DSP.unwrap(angle.(field.ψ[end,:]),range=pi)
φₜₕ = γ*abs2.(ψ₀)*L

@test sum(abs.((φₜₕ.-φ)./φₜₕ))/length(φₜₕ) < 0.001

end
