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
include("raman.jl")
include("simulation.jl")
include("analysis.jl")

end
