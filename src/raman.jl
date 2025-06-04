"""
Agrawaal Raman model

Default values are (fr = 0.18 , τl = 32e-15 [s], τvib = 12.2e-15 [s])

Construct the Agrawal Raman model as a [`RamanModel`](@ref) with typical silica fiber parameters.

# Arguments
- `fr`: Fractional Raman contribution (default: 0.18)
- `τl`: Raman longitudinal relaxation time in seconds (default: 32e-15)
- `τvib`: Raman vibrational period in seconds (default: 12.2e-15)

# Returns
- `RamanModel`: A RamanModel struct with the Agrawal time response.

# References
- G. P. Agrawal, Nonlinear Fiber Optics, 5th Edition, Eq. (7.1.3)

# See also
- [`RamanModel`](@ref)
- [`NoRaman`](@ref)
"""
function raman_linagrawaal(fr = 0.18, τl = 32e-15,
	τvib = 12.2e-15)

	function time_response(t)
		hr = 0.0 .+ (t .>= 0) .* ((τvib^2 + τl^2) / τvib / τl^2 * exp.(-t / τl) .* sin.(t / τvib))
		return hr
	end

	return RamanModel(fr, time_response)

end
