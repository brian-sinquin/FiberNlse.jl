module FiberNlse
using FFTW

# Physical units & constants

nm = ns = 1e-9
ps = pm = 1e-12
km = kHz = 1e3
mW = mm = ms = 1e-3
GHz = 1e9
THz = 1e12
m = W = 1


# light speed in vaccum
c = 299792458 # m/s

export Fiber, Dispersion, Field, dispersion, smf28, edfa, propagate, output, concatf
export ps,pm,ns,nm,mW,mm,ms,m,W,km,kHz,GHz,THz

include("spacetime.jl")
include("dispersion.jl")
include("fiber.jl")
include("raman.jl")
include("simulation.jl")
include("analysis.jl")

end
