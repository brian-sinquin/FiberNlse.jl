module FiberNlse

using FFTW

# Export only the main public API
export Waveguide, Solution, GNLSEProblem, NoRaman, RamanModel
export gnlse, combine, input, output, raman_linagrawaal

# Light celerity
c = 299792458

@info "The integrated signal is supposed to follow the negative phase convention."

include("datatypes.jl")
include("api.jl")
include("raman.jl")
include("nonlinearity.jl")
include("utils.jl")
include("solvers.jl")
end

