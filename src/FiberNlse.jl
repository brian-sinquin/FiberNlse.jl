module FiberNlse
using DSP: fft, ifft, ifftshift, fftshift, unwrap, rms

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

export Fiber,
    Dispersion, Field, dispersion, smf28, edfa, propagate, output, concatf, phase, instFreq
export ps, pm, ns, nm, mW, mm, ms, m, W, km, kHz, GHz, THz

global DefaultPhaseConvention::Symbol
global DefaultPhaseConvention = :positive

function setPhaseConvention(conv::Symbol)
    global DefaultPhaseConvention = conv
end

getPhaseConvention() = DefaultPhaseConvention
function getPhaseConventionSign(conv::Symbol)
    conv == :positive ? +1 : -1
end


include("spacetime.jl")
include("dispersion.jl")
include("fiber.jl")
include("simulation.jl")
include("analysis.jl")
include("utils.jl")

end
