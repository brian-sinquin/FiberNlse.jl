include("../src/FiberNlse.jl")
using Plots;gr()
using FFTW;
using Peaks;
using SpecialFunctions;
U = FiberNlse; # Alias

begin
    
# CONSTANTS
    T = 50000*U.ps;
    Fs = 200*U.GHz;
# Simulation dimension
    Nₜ, Nₗ = (Int(round(Fs*T)),500);

# Fiber properties

    D = -7*U.ps/U.nm/U.km;
    α = 0*0.046/U.km;
    γ =  1.1/U.W/U.km;
    L = 4131;
    fib = U.Fiber(L, α, D, γ)


# Signal properties
λ1 = 1550.0*U.nm; # Wavelength
λ2 = λ1+1*U.pm; # Wavelength
λ0 = 0.5*(λ1+λ2);
Δλ = 0.5*abs(λ1-λ2)
sim, t = U.configure(Nₜ,Nₗ,D, γ, α, L, T, λ0);
Δw = 2pi*U.c*abs(λ1-λ2)/(λ1*λ2)
end

SP(ψ) = fftshift(abs.(fft(ψ)).^2)/sim.Nₜ^2 
ν = fftshift(FFTW.fftfreq(sim.Nₜ, 1.0 / (sim.dt)));

function run(P)
    # Input construction
    Ψₒ = @. sqrt(0.5*P)*(exp.(1im*Δw.*t).+exp.(-1im*Δw.*t))
    U.inputSignal(sim,Ψₒ);
    U.simulate(sim, false);
    return (sim.Ψ[1,:],sim.Ψ[end,:])
end


ψi,ψf = run(0.5);
sp = SP(ψf)
idx = Peaks.argmaxima(sp)
idx = idx[findall(ν[idx] .< 4*Δw/(2pi))]
idx = idx[findall(ν[idx] .> 0)]
plot(SP(ψf))
scatter!(idx, sp[idx])

begin
Pl = (30:10:300)*U.mW
r = []
pf = []
for p in Pl
ψi,ψf = run(p);
sp = SP(ψf)
idx = Peaks.argmaxima(sp)
idx = idx[findall(ν[idx] .< 4*Δw/(2pi))]
idx = idx[findall(ν[idx] .> 0)]

if length(idx)>1
    append!(r, 10*log10(sp[idx[1]]/sp[idx[2]]))
    append!(pf, p)
end
end
end
pth = range(30, 300, 1000)*U.mW
a=0.82pi/maximum(pth)

rr(x) = @. (besselj0(0.5x)^2+besselj1(0.5x)^2)/(besselj1(0.5x)^2+besseljx(2,0.5x)^2)
scatter(pf,r)
plot!(pth, 10*log10.(rr.((a*pth))))


γexp = 0.5*a/L
δγ = abs(γexp-γ)/γ*100