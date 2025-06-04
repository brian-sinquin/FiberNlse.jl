"""
    _erk4ip_step!(stepper, model)

Perform a single embedded Runge-Kutta 4(5) IP step for the GNLSE integration.

This function updates the fields in `stepper` in-place using the current model parameters, including the nonlinear operator and dispersion term.

# Arguments
- `stepper::Stepper`: The stepper struct holding the current state and intermediate arrays.
- `model::GNLSEProblem`: The problem definition containing all precomputed arrays and parameters.

# Modifies
- Updates all relevant fields of `stepper` in-place for the next integration step.

# References
- Hult, J. (2007). "A Fourth-Order Runge–Kutta in the Interaction Picture Method for Simulating Supercontinuum Generation in Optical Fibers." J. Lightwave Technol. 25, 3770-3775. [doi:10.1109/JLT.2007.909373](https://doi.org/10.1109/JLT.2007.909373)

# See also
- [`Stepper`](@ref)
- [`GNLSEProblem`](@ref)
"""
function _erk4ip_step!(stepper, model)

	#todo! inplace rk_order + nl fct

	stepper.e = exp.(0.5 * stepper.dz * model.dispersion_term)

	stepper.Uip = stepper.e .* stepper.U

	stepper.k1 = stepper.e .* stepper.NU

	stepper.k2 = model.nonlinear_function(model.ifftp * (stepper.Uip .+ 0.5 * stepper.dz * stepper.k1), model)

	stepper.k3 = model.nonlinear_function(model.ifftp * (stepper.Uip .+ 0.5 * stepper.dz * stepper.k2), model)

	stepper.k4 = model.nonlinear_function(model.ifftp * (stepper.e .* (stepper.Uip .+ stepper.dz * stepper.k3)), model)

	stepper.r = stepper.e .* (stepper.Uip .+ stepper.dz * (stepper.k1 / 6.0 .+ stepper.k2 / 3.0 .+ stepper.k3 / 3.0))

	stepper.U1 = stepper.r .+ stepper.dz * stepper.k4 / 6.0

	stepper.k5 = model.nonlinear_function(model.ifftp * stepper.U1, model)

	stepper.U2 = stepper.r .+ stepper.dz * (stepper.k4 / 15.0 .+ stepper.k5 / 10.0)

end
