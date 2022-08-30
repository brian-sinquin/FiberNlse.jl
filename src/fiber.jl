"""
    Fiber(L::Float64, D::Dispersion, γ::Float64, α::Float64)

    Structure containing fiber caracteristics :
    - L : length in [m]
    - D : (dispersion)[@ref Dispersion]
    - γ : nonlinear factor in [1 / W / m]
    - α : fiber attenuation or gain (exponential) in [1 / m]

    See also : (smf28)[@ref smf28], (edfa)[@ref edfa]
"""
struct Fiber
    L::Float64
    D::Dispersion
    γ::Float64
    α::Float64
    λ::Float64
end

"""
    edfa(L::Float64, G::Float64, λ::Float64)

    Generates a fiber of length L [m] and gain G [dB] arround a wavelength of λ [m]

    See also : (smf28)[@ref smf28]
"""
function edfa(L::Float64, G::Float64, λ::Float64)
    g  = 0.1*log(10)*G/L
    Fiber(L,disperion(-10e-6,λ), 1.1e-3, -g, λ);
end

"""
    smf28(L::Float64, λ::Float64)

    Generates a standard monomode fiber (smf28) of length L [m] arround a wavelength of λ [m]

    See also : (edfa)[@ref edfa]
"""
function smf28(L::Float64, λ::Float64)
    return Fiber(L, dispersion(17e-6, λ), 1.1e-3, 0.046e-3, λ)
end
