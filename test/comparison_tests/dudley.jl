using Test
using CSV
using DataFrames
using FiberNlse

@testset "Supercontinuum Generation - Dudley Test" begin
    # === Simulation Parameters ===
    n = 2^13                   # Number of grid points
    twidth = 12.5e-12          # Width of time window [s]
    c = 299792458              # Speed of light [m/s]
    wavelength = 835e-9        # Reference wavelength [m]
    w0 = (2 * pi * c) / wavelength   # Reference frequency [Hz]
    dt = twidth / n
    T = (-n/2:n/2-1) .* dt     # Time grid

    # === Input Pulse ===
    power = 10000              # Peak power of input [W]
    t0 = 28.4e-15              # Duration of input [s]
    A = @. sqrt(power) * sech(T / t0) .+ 0.0im  # Input field [W^(1/2)]

    # === Fiber Parameters ===
    flength = 0.15             # Fiber length [m]
    betas = [0.0, -1.1830e-026, 8.1038e-041, -9.5205e-056, 2.0737e-070,
        -5.3943e-085, 1.3486e-099, -2.5495e-114, 3.0524e-129,
        -1.7140e-144]     # Dispersion coefficients [s^n/m]
    gamma = 0.11               # Nonlinear coefficient [1/W/m]
    loss = 0.0                 # Loss [dB/m]

    # Define the waveguide
    wg = Waveguide(
        loss,
        betas,
        gamma,
        wavelength,
        flength,
        raman_model=raman_linagrawaal(),
        self_steepening=true
    )

    # === Problem Setup ===
    nsaves = 200               # Number of length steps to save field at
    dz = flength / (0.5 * nsaves)  # Step size for propagation

    # Solve the GNLSE
    sol = gnlse(A, T, wg, nsaves=nsaves, dz=dz, reltol=1e-6)

    # === Load Reference Data ===
    fn = joinpath(dirname(@__FILE__), "data/table_dudley_test_t.csv")
    dat = CSV.read(fn, DataFrame)

    t_dudley = dat.t
    It_dudley = parse.(ComplexF64, dat.At) .|> abs2

    # === Compare Results ===
    I = reverse(abs2.(sol.At[end, :]))  # Reverse the simulated intensity
    err = 1 / length(I) * sum(abs.(I .- It_dudley) / maximum(I))  # Compute error

    # Test if the error is within acceptable bounds
    @test err < 1.0 / 100
end
