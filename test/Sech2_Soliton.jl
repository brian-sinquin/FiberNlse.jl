include("../src/FiberNlse.jl")
using .FiberNlse

# Simulation dimension
Nₜ, Nₗ = (1000,1000);

# Fiber properties
L = 5.0*km; # Fiber length
fib = FiberNlse.smf28(L)
# Signal properties
T = 40*ps; # Signal duration
λ = 1550*nm; # Wavelength
τ = 5*ps; # Pulse duration
N = 1 # Soliton number
sim,t,l = FiberNlse.configure(Nₜ,Nₗ,fib, T, λ);

# Input construction
Pp = abs((sim.β2/fib.γ/τ^2)*N^2) # Soliton power
Ψₒ = @. sqrt(Pp)/cosh(t/τ) # Soliton formula
FiberNlse.inputSignal(sim,Ψₒ);

FiberNlse.simulate(sim, false); # run the simulation

# Visualization

using Plots
heatmap(t/ps, l/km, abs.(sim.Ψ').^2)

#= using GLMakie

begin
fig = Figure(resolution=(860, 400), fontsize=12)
ax1 = Axis3(fig[1,1]; aspect=(1,1,1),ylabel="Distance [km]", xlabel="Time [ps]", zlabel="Power [W]")
ax2 = Axis(fig[1,2], ylabel="Distance [km]", xlabel="Time [ps]")
hm = surface!(ax1, t/ps, l/km, abs.(sim.Ψ').^2)
heatmap!(ax2, t/ps, l/km, abs.(sim.Ψ').^2)
Colorbar(fig[1, 3], hm, height=Relative(0.5), label="Power [W]")
fig
end
 =#