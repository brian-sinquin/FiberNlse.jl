#TODO!: inplace functions

"""
    _spm(u, model)

Nonlinear operator for self-phase modulation (SPM) only (no Raman, no self-steepening).

# Arguments
- `u`: Input field (time domain).
- `model`: GNLSEProblem or model struct with parameters.

# Returns
- Nonlinear term in frequency domain.
"""
function _spm(u, model)
	return 1.0im .* model.γ * (model.fftp * (u .* abs2.(u)))
end

"""
    _spm_self_steepening(u, model)

Nonlinear operator for SPM with self-steepening (no Raman).

# Arguments
- `u`: Input field (time domain).
- `model`: GNLSEProblem or model struct with parameters.

# Returns
- Nonlinear term in frequency domain.
"""
function _spm_self_steepening(u, model)
	return 1.0im .* model.γ * (1 .+ model.ω / model.ω0) .* (model.fftp * (u .* abs2.(u)))
end

"""
    _spm_self_steepening_raman(u, model)

Nonlinear operator for SPM with self-steepening and Raman effect.

# Arguments
- `u`: Input field (time domain).
- `model`: GNLSEProblem or model struct with parameters.

# Returns
- Nonlinear term in frequency domain.
"""
function _spm_self_steepening_raman(u, model)
	IT = abs2.(u)
	RS = model.dt * model.fr * (model.ifftp * ((model.fftp * IT) .* model.raman_freq_response))
	return 1.0im .* model.γ * (1 .+ model.ω / model.ω0) .* (model.fftp * (u .* ((1.0 - model.fr) .* IT .+ RS)))
end

"""
    _spm_raman(u, model)

Nonlinear operator for SPM with Raman effect (no self-steepening).

# Arguments
- `u`: Input field (time domain).
- `model`: GNLSEProblem or model struct with parameters.

# Returns
- Nonlinear term in frequency domain.
"""
function _spm_raman(u, model)
	IT = abs2.(u)
	RS = model.dt * model.fr * (model.ifftp * ((model.fftp * IT) .* model.raman_freq_response))
	return 1.0im .* model.γ .* (model.fftp * (u .* ((1.0 - model.fr) .* IT .+ RS)))
end

"""
    choose_nonlinear_term(self_steepening::Bool = true, raman::Bool = true)

Select the appropriate nonlinear operator function based on self-steepening and Raman effect flags.

# Arguments
- `self_steepening::Bool`: Whether to include self-steepening.
- `raman::Bool`: Whether to include Raman effect.

# Returns
- Function: The selected nonlinear operator function.

# See also
- [`_spm`](@ref)
- [`_spm_self_steepening`](@ref)
- [`_spm_self_steepening_raman`](@ref)
- [`_spm_raman`](@ref)
"""
function choose_nonlinear_term(self_steepening::Bool = true, raman::Bool = true)

	if self_steepening && raman
		return _spm_self_steepening_raman
	elseif self_steepening && !raman
		return _spm_self_steepening
	elseif !self_steepening && raman
		return _spm_raman
	else
		return _spm
	end

end
