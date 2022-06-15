include("../src/FiberNlse.jl")
using .FiberNlse
# Simulation dimension
Nₜ, Nₗ = (1000,1000);

# Fiber properties
L = 5km; # Fiber length
fib2 = FiberNlse.smf28(2*L)
fib1 = FiberNlse.Fiber(-fib1.D, fib1.α, fib1.γ,L)

# Signal properties
T = 200*ps; # Signal duration
λ = 1550*nm; # Wavelength
τ = 5*ps; # Pulse duration
N = 1 # Soliton number
sim,t,l = FiberNlse.configure(Nₜ,Nₗ,fib1, T, λ);
sim2,t,l = FiberNlse.configure(Nₜ,Nₗ,fib2, T, λ);
# Input construction
Pp = abs((sim.β2/fib1.γ/τ^2)*N^2) # Soliton power
Ψₒ = @. sqrt(Pp)/cosh(t/τ) # Soliton formula
FiberNlse.inputSignal(sim,Ψₒ);
FiberNlse.simulate(sim, true); # run the simulation
FiberNlse.transition(sim,sim2)
FiberNlse.simulate(sim2, true);
# Visualization
S = vcat(sim.Ψ,sim2.Ψ)
using Plots
gr()
heatmap(abs2.(S'))

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