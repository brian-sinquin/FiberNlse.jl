module FiberNlse
using FFTW

# Physical units & constants

nm = ns = 1e-9
ps = pm = 1e-12
km = kHz = 1e3
mW = mm = 1e-3
GHz = 1e9
Thz = 1e12
m = 1
W = 1

# light speed in vaccum
c = 299792458

export Fiber, Dispersion, Field, dispersion, smf28, edfa, propagate


include("spacetime.jl")
include("dispersion.jl")
include("fiber.jl")
include("simulation.jl")
include("analysis.jl")


# chembo
#= Nₗ=1000
λ = 1550e-9
t = LinRange(0,10000, 10000)*1e-12
f0 = 1e10
P₀=0.5
Vpidc=6.6
Vpirf=4.2
V0=Vpirf/2
Vpip=4.
θ=pi
Vb=Vpidc/2
q=1.6
y = @. sqrt(P₀)*cos(0.5*pi*(V0/Vpirf*cos(2pi*f0*t)+Vb/Vpidc))*exp(1im*pi*q*V0/Vpip*cos(2pi*f0*t+θ))
fib = Fiber(4km,dispersion(20*ps/nm/km, λ),1.1 / W/km, 0.046/km)
field = simulate(y, fib,t[end], Nₗ)
using DSP
heatmap(abs2.(field.ψ))
plot(abs2.(field.ψ[end,:]))
plot(10log10.(abs2.(fftshift(fft(field.ψ[end,:])))))
plot(unwrap(angle.(fftshift(fft(field.ψ[end,:])))))
xlims!(300,700) =#
end