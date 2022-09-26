include("FiberNlse.jl")
using .FiberNlse
using Plots, DSP

ps = 1e-12
t = (-100:0.1:100) * ps
T = abs(t[1] - t[end])
t0 = 2ps

E = @. 1e-10 * exp(-(t / t0)^2)

f0 = 500e9
L = 1000.0
λ = 1550e-9
fib = FiberNlse.smf28(L, λ)
fib = Fiber(L, Dispersion(-2.16e-26), 0, 0, λ)

field = FiberNlse.propagate(E, fib, T, 500)

Ef = output(field)

plot(t, abs.(E))
plot!(t, abs.(Ef))

function deriv(y::AbstractVector, x::AbstractVector)
    function centraldiff(v::AbstractVector)
        dv = diff(v) / 2 # half the derivative
        a = [dv[1]; dv] # copies first element
        a .+= [dv; dv[end]] # copies last element, add both results to compute average
        return (a)
    end
    return centraldiff(y) ./ centraldiff(x)
end


phi = -deriv(unwrap(angle.(E)), t) / 2pi / 1e9
phi2 = -deriv(unwrap(angle.(Ef)), t) / 2pi / 1e9

plot(t, phi)
plot!(t, phi2)

plot(t, abs2.((E .+ 1 * exp.(im * 2pi * f0 * t))))
plot!(t, abs2.(Ef .+ 1 * exp.(im * 2pi * f0 * t)))
