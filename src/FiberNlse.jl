module FiberNlse
using DSP: fft, ifft, ifftshift, fftshift, unwrap, rms
using DifferentialEquations

# Physical units & constants

nm = ns = 1e-9
ps = pm = 1e-12
km = kHz = 1e3
mW = mm = ms = 1e-3
GHz = 1e9
THz = 1e12
m = W = Hz = 1

# light speed in vaccum
c = 299792458 # m/s


global phase_convention::Symbol = :positive




export Fiber,
    Dispersion, Field, dispersion, smf28, edfa, propagate, output, concatf, phase, instFreq, propagate4,derivate
export ps, pm, ns, nm, mW, mm, ms, m, W, km, kHz, GHz, THz


function setPhaseConvention(conv::Symbol)
    global phase_convention = conv
end

getPhaseConvention() = phase_convention
getPhaseConventionSign() = phase_convention == :positive ? +1 : -1
getPhaseConventionSign(symbol::Symbol) = symbol == :positive ? +1 : -1

include("core.jl")
include("spacetime.jl")
include("dispersion.jl")
include("fiber.jl")
include("simulation.jl")
include("analysis.jl")
include("utils.jl")

end
