"""
	RamanModel

Struct describing the Raman response of the waveguide.

# Fields
- `fr::Float64`: Fractional contribution of the Raman effect to the overall nonlinear effect.
- `time_response::Function`: User-defined function returning the Raman impulse response as a function of time.

# See also
- [`raman_linagrawaal`](@ref)
- [`NoRaman`](@ref)
- [`Waveguide`](@ref)
"""
struct RamanModel
	fr::Float64
	time_response::Function
end

"""
	NoRaman

An empty [`RamanModel`](@ref) struct to avoid simulating the Raman effect.

# See also
- [`RamanModel`](@ref)
- [`Waveguide`](@ref)
"""
const NoRaman = RamanModel(0.0, () -> ())



"""
	Waveguide

Struct describing the waveguide, i.e., the propagation conditions.

# Fields
- `α::Union{Float64, Array{Float64}}`: Linear optical intensity loss factor (scalar or frequency dependent)
- `βs::Array{Float64}`: GVD orders vector (starting at 2nd order)
- `γ::Union{Float64, Array{Float64}}`: Nonlinear factor (scalar or frequency dependent)
- `λc::Float64`: Center wavelength of the problem (signal and dispersion)
- `L::Float64`: Physical length of the waveguide
- `raman_model::RamanModel`: Raman model of the waveguide ([`RamanModel`](@ref))
- `self_steepening::Bool`: Whether to take into account the self-steepening effect

# See also
- [`RamanModel`](@ref)
- [`NoRaman`](@ref)
- [`GNLSEProblem`](@ref)
"""
mutable struct Waveguide
	α::Union{Float64, Array{Float64}}
	βs::Array{Float64}
	γ::Union{Float64, Array{Float64}}
	λc::Float64
	L::Float64
	raman_model::RamanModel
	self_steepening::Bool
end

"""
	Waveguide(α, βs, γ, λc, L; raman_model=NoRaman, self_steepening=false)

Constructs a [`Waveguide`](@ref) struct.

By default, no Raman effect nor self-steepening is activated.

# See also
- [`Waveguide`](@ref)
- [`RamanModel`](@ref)
- [`NoRaman`](@ref)
"""
Waveguide(α, βs, γ, λc, L; raman_model = NoRaman, self_steepening = false) = Waveguide(α, βs, γ, λc, L, raman_model, self_steepening)

"""
	GNLSEProblem

Struct containing all the parameters and precomputed arrays needed for solving the generalized nonlinear Schrödinger equation (GNLSE).

# Fields
- `ω::AbstractArray{Float64}`: Frequency grid
- `dt::Float64`: Time step
- `N::Int`: Number of time/frequency points
- `fftp::Any`: FFT plan
- `ifftp::Any`: IFFT plan
- `dispersion_term::AbstractArray{ComplexF64}`: Precomputed dispersion term
- `nonlinear_function::Function`: Nonlinear operator function
- `fr::Float64`: Fractional Raman contribution
- `γ::Float64`: Nonlinear coefficient
- `raman_freq_response::Union{Nothing, AbstractArray{ComplexF64}}`: Raman frequency response (if any)
- `ω0::Float64`: Central frequency
- `L::Float64`: Waveguide length

# See also
- [`Waveguide`](@ref)
- [`gnlse`](@ref)
- [`Solution`](@ref)
"""
mutable struct GNLSEProblem
	ω::AbstractArray{Float64}
	dt::Float64
	N::Int
	fftp::Any
	ifftp::Any
	dispersion_term::AbstractArray{ComplexF64}
	nonlinear_function::Function
	fr::Float64
	γ::Float64
	raman_freq_response::Union{Nothing, AbstractArray{ComplexF64}}
	ω0::Float64
	L::Float64
end


"""
	Solution

Struct containing the result of the simulation.

# Fields
- `z::AbstractArray{Float64}`: Distance vector
- `t::AbstractArray{Float64}`: Time vector
- `f::AbstractArray{Float64}`: Frequency vector
- `At::Matrix{ComplexF64}`: Time-domain propagation field
- `Af::Matrix{ComplexF64}`: Frequency-domain (spectrum) propagation field

# See also
- [`combine`](@ref)
- [`GNLSEProblem`](@ref)
"""
mutable struct Solution
	z::AbstractArray{Float64}
	t::AbstractArray{Float64}
	f::AbstractArray{Float64}
	At::Matrix{ComplexF64}
	Af::Matrix{ComplexF64}
end

"""
	Stepper

Struct containing the fields that are updated or needed during the solver iterations.

# Fields
- `U::AbstractArray{ComplexF64}`: Current field (frequency domain)
- `NU::AbstractArray{ComplexF64}`: Nonlinear operator applied to field
- `dz::Float64`: Current step size
- `z::Float64`: Current propagation distance
- `local_error::Float64`: Local error estimate
- `k1`, `k2`, `k3`, `k4`, `k5`: Runge-Kutta intermediate steps
- `Uip`, `U1`, `U2`: Intermediate fields
- `e`, `r`: Temporary arrays for computation
- `it::Int`: Iteration counter

# See also
- [`GNLSEProblem`](@ref)
- [`gnlse`](@ref)
"""
mutable struct Stepper
	U::AbstractArray{ComplexF64}
	NU::AbstractArray{ComplexF64}
	dz::Float64
	z::Float64
	local_error::Float64
	k1::AbstractArray{ComplexF64}
	k2::AbstractArray{ComplexF64}
	k3::AbstractArray{ComplexF64}
	k4::AbstractArray{ComplexF64}
	k5::AbstractArray{ComplexF64}
	Uip::AbstractArray{ComplexF64}
	U1::AbstractArray{ComplexF64}
	U2::AbstractArray{ComplexF64}
	e::AbstractArray{ComplexF64}
	r::AbstractArray{ComplexF64}
	it::Int
end


