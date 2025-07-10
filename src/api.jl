"""
    GNLSEProblem(t::AbstractArray{Float64}, wg::Waveguide)

Constructs a `GNLSEProblem` struct containing all parameters and precomputed arrays needed for solving the generalized nonlinear Schrödinger equation (GNLSE) for a given time grid `t` and `Waveguide` `wg`.

# See also
- [`Waveguide`](@ref)
- [`gnlse`](@ref)
"""
function GNLSEProblem(t::AbstractArray{Float64}, wg::Waveguide)

	fftp = plan_fft(t, flags = FFTW.MEASURE)
	ifftp = plan_ifft(t, flags = FFTW.MEASURE)

	T = t[end] - t[1]
	N = length(t)

	ω = 2π .* fftshift((-N÷2:N÷2-1) ./ T)

	# raman model data / no value if no raman 
	raman_freq_response = nothing
	if wg.raman_model.fr != 0.0
		raman_freq_response = conj((fftp * ifftshift(wg.raman_model.time_response(t))))
	end

	nonlinear_function = choose_nonlinear_term(wg.self_steepening, !isnothing(raman_freq_response))
	dispersion_term = -0.5wg.α .+ 1im * sum([(wg.βs[i] / factorial(i)) .* ω .^ i for i in eachindex(wg.βs)])

	GNLSEProblem(ω, T / N, N, fftp, ifftp, dispersion_term, nonlinear_function, wg.raman_model.fr, wg.γ, raman_freq_response, 2pi * c / wg.λc, wg.L)
end


"""
    _compute_error(a, b)

Computes the normalized root mean square error between two arrays `a` and `b`.

# See also
- [`_compute_error!`](@ref)
"""
function _compute_error(a, b)
	return sqrt(sum(abs2.(a .- b)) ./ sum(abs2.(a)))
end


"""
    _compute_error!(stepper::Stepper, a, b)

Updates the `local_error` field of `stepper` with the normalized root mean square error between arrays `a` and `b`.

# See also
- [`_compute_error`](@ref)
"""
function _compute_error!(stepper::Stepper, a, b)
	stepper.local_error = sqrt(sum(abs2.(a .- b)) ./ sum(abs2.(a)))
end

"""
    _integrate_to_z(stepper::Stepper, z::Float64, prob::GNLSEProblem, maxiters::Int, reltol::Float64)

Integrates the GNLSE from the current position of `stepper` up to distance `z`, using adaptive step size control. Throws an error if `maxiters` is exceeded.

# See also
- [`gnlse`](@ref)
"""
function _integrate_to_z(stepper::Stepper, z::Float64, prob::GNLSEProblem, maxiters::Int, reltol::Float64)
	stepper.it = 0
	while stepper.z < z
		_erk4ip_step!(stepper, prob)
		_compute_error!(stepper, stepper.U1, stepper.U2)

		dzopt =
			max(0.5, min(2.0, 0.9 * sqrt(sqrt(reltol / stepper.local_error)))) * stepper.dz

		if stepper.local_error <= reltol
			stepper.dz = min(dzopt, abs(z - stepper.z))
			stepper.z += stepper.dz
			stepper.U = stepper.U1
			stepper.NU = stepper.k5
		else
			stepper.dz = dzopt
			stepper.it += 1
			if (stepper.it >= maxiters)
				throw(ErrorException("Max number of iteration exceeded!"))
			end
		end
	end
end

"""
    gnlse(u::AbstractArray{ComplexF64}, t::AbstractArray{Float64}, prob::GNLSEProblem; nsaves=20, dz=1.0, reltol=1e-6, maxiters=1000)

Solves the GNLSE for initial field `u` and time grid `t` using the problem definition `prob`. Returns a `Solution` struct with the propagation results.

# See also
- [`GNLSEProblem`](@ref)
- [`Solution`](@ref)
"""
function gnlse(u::AbstractArray{ComplexF64}, t::AbstractArray{Float64}, prob::GNLSEProblem; nsaves = 20, dz = 1.0, reltol = 1e-6, maxiters = 1000)

	dz = min(prob.L / (2 * nsaves), dz)
	k_init = similar(u)
	# initial stepper
	stepper = Stepper(prob.fftp * u, prob.nonlinear_function(u, prob), dz, 0.0, 0.0, k_init, k_init, k_init, k_init, k_init, k_init, k_init, k_init, k_init, k_init, 0)
	zsaves = (0:nsaves) * prob.L / nsaves
	M = zeros(ComplexF64, (nsaves + 1, prob.N))
	M[1, :] = stepper.U
	ϵ_hist = zeros(nsaves + 1)


	for i ∈ 2:nsaves+1
		_integrate_to_z(stepper, zsaves[i], prob, maxiters, reltol)
		ϵ_hist[i] = stepper.local_error
		M[i, :] = stepper.U
	end

	return Solution(zsaves, t, prob.ω / 2pi, ifft(M, 2), M)
end

"""
    gnlse(u::AbstractArray{ComplexF64}, t::AbstractArray{Float64}, wg::Waveguide; args...)

Convenience wrapper for `gnlse` that constructs a `GNLSEProblem` from a `Waveguide` and calls `gnlse`.

# See also
- [`GNLSEProblem`](@ref)
"""
gnlse(u::AbstractArray{ComplexF64}, t::AbstractArray{Float64}, wg::Waveguide; args...) = gnlse(u, t, GNLSEProblem(t, wg); args...)

"""
    gnlse(u::AbstractArray{ComplexF64}, t::AbstractArray{Float64}, probs::Vector{GNLSEProblem}; args...)

Propagates the field `u` through a sequence of `GNLSEProblem`s, chaining the output of each as the input to the next. Returns a combined `Solution`.

# See also
- [`combine`](@ref)
"""
function gnlse(u::AbstractArray{ComplexF64}, t::AbstractArray{Float64}, probs::Vector{GNLSEProblem}; args...)
	sols = [gnlse(u, t, probs[1]; args...)]
	for prob ∈ probs[2:end]
		push!(sols, gnlse(sols[end].At[end, :], t, prob; args...))
	end
	combine(sols)
end

"""
    gnlse(u::AbstractArray{ComplexF64}, t::AbstractArray{Float64}, wgs::Vector{Waveguide}; args...)

Propagates the field `u` through a sequence of `Waveguide`s, constructing a `GNLSEProblem` for each. Returns a combined `Solution`.

# See also
- [`GNLSEProblem`](@ref)
"""
gnlse(u::AbstractArray{ComplexF64}, t::AbstractArray{Float64}, wgs::Vector{Waveguide}; args...) = gnlse(u, t, [GNLSEProblem(t, wg) for wg in wgs]; args...)

