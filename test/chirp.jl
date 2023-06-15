using FiberNlse, Plots, StatsBase, DSP
FiberNlse.setPhaseConvention(:negative)
D = 17e-6
α = 0.026e-3
γ = 0*1.1e-3
L = 5000.0

# Signal properties
λ = 1550e-9 # Wavelength

Nt, Nz = (2^10,200)
T = 300e-12
t0 = 30e-12
t = LinRange(-1,1,Nt)*T/2

LEF = -3.8
P0 =0.08
y = sqrt(P0).*exp.(-((t./t0).^2))
y= y.*cis.(0.5*LEF*log.(abs2.(y)./abs2.(y[1])))
S = sqrt(P0)*y./sqrt(mean(abs2.(y)))


fib = Fiber(L, dispersion(-D, λ), γ, α, λ)
field = propagate(S, fib, T, Nt)

φ = instFreq(field.ψ[end, :],t)/1e9
φₜₕ = instFreq(S, t)/1e9

begin
    plot(t/ps,φ, label="end")
    plot!(t/ps,φₜₕ, label="begin")

end

begin
    plot(t/ps,field.ψ[end, :] .|> abs2, label="end")
    plot!(t/ps,S .|> abs2, label="begin")
    xlims!(-50,50)
end

heatmap(abs2.(field.ψ))