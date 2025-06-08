# # Supercontinuum Generation Example
#
# This example demonstrates supercontinuum generation in an optical fiber using the GNLSE.

# Import necessary packages
using FiberNlse
using Plots
using FFTW

# Simulation parameters
n = 2^13                   # number of grid points
twidth = 12.5e-12          # width of time window [s]
c = 299792458              # speed of light [m/s]
wavelength = 835e-9        # reference wavelength [m]
w0 = (2 * pi * c) / wavelength   # reference frequency [Hz]
dt = twidth / n
T = (-n/2:n/2-1) .* dt # time grid

# Input pulse
power = 10000              # peak power of input [W]
t0 = 28.4e-15              # duration of input [s]
A = @. sqrt(power) * sech(T / t0) .+ 0.0im # input field [W^(1/2)]

# Fiber parameters
flength = 0.15             # fibre length [m]
betas = [0.0, -1.1830e-026, 8.1038e-041, -9.5205e-056, 2.0737e-070,
    -5.3943e-085, 1.3486e-099, -2.5495e-114, 3.0524e-129,
    -1.7140e-144]
gamma = 0.11               # nonlinear coefficient [1/W/m]
loss = 0.0                   # loss [dB/m]

# Define the waveguide
wg = Waveguide(loss=loss, βs=betas, γ=gamma, raman_model=raman_linagrawaal(), self_steepening=true)

# Simulation setup
nsaves = 200     # number of length steps to save field at
solver = Solver(flength, nsaves)

# Solve the GNLSE
sol = gnlse(A, wg, T, solver, reltol=1e-6)

# Plotting the results
z = sol.z
It = abs2.(sol.A)

plt_t = plot(T * 1e12, It[1, :], legend=false,
    xlabel="Time (ps)", ylabel="Intensity (W)",
    title="Supercontinuum Generation")

gif = @animate for i in 1:length(z)
    plot(T * 1e12, It[i, :], legend=false,
        xlabel="Time (ps)", ylabel="Intensity (W)",
        title="Supercontinuum Generation",
        xlims=(-5, 5), ylims=(0, 1.1 * maximum(It)))
end

gifpath = joinpath(dirname(@__FILE__), "supercontinuum.gif")
gif(gif, gifpath, fps=10)

println("Supercontinuum simulation complete.  Check the folder ", dirname(@__FILE__), " for the supercontinuum.gif")

# Display the final time slice
display(plt_t)
