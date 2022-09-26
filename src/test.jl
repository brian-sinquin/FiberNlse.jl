using FiberNlse
using Plots, DSP

ps = 1e-12
t = (-100:0.1:100) * ps
T = (t[1] - t[end]) |> abs 
t0 = 20ps
P0 = 0.5
E = @. sqrt(P0) * exp(-(t / t0)^2)
D = 18e-6
f0 = 500e9
L = 1000.0
λ = 1550e-9
γ = 1.1e-3
fib = Fiber(L, dispersion(D, λ), γ, 0, λ)
field = FiberNlse.propagate(E, fib, T, 500)

Ef = output(field)

plot(t, abs2.(E))
plot!(t, abs2.(Ef))
begin
    finsti = instFreq(E,t)/1e9
    finstf = instFreq(Ef,t)/1e9
    plot(t, finsti, label="begin")
    plot!(t, finstf, label="end")
end