# # Soliton Propagation Example
#
# This example demonstrates the propagation of a fundamental soliton in an optical fiber.

# Import necessary packages
using FiberNlse
using Plots

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
fib = Waveguide(α=α, β2=-2.6e-26, γ=1.1e-3)

t = range(-T / 2, T / 2, length=N)

# Input construction
P₀ = abs((fib.β2 / fib.γ / τ^2) * n^2) # Soliton power
Ψₒ = sqrt(P₀) * sech.(t ./ τ) .+ 0.0im # Soliton formula

# Create a Solver object
solver = Solver(L, 200)

# Propagate the pulse
sol = gnlse(Ψₒ, fib, t, solver)

# Plotting the results
z = sol.z
It = abs2.(sol.A)

plt_t = plot(t * 1e12, It[1, :], legend=false,
    xlabel="Time (ps)", ylabel="Intensity (W)",
    title="Soliton Propagation")

gif = @animate for i in 1:length(z)
    plot(t * 1e12, It[i, :], legend=false,
        xlabel="Time (ps)", ylabel="Intensity (W)",
        title="Soliton Propagation",
        xlims=(-5 * τ * 1e12, 5 * τ * 1e12), ylims=(0, 1.1 * maximum(It)))
end

gifpath = joinpath(dirname(@__FILE__), "soliton.gif")
gif(gif, gifpath, fps=10)

println("Soliton simulation complete. Check the folder ", dirname(@__FILE__), " for the soliton.gif")

# Display the final time slice
display(plt_t)
