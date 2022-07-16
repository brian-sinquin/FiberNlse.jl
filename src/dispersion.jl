"""
    Dispersion(βₛ::Union{Float64,Vector{Float64}})

    Structure containing dispersion caracteristics for the fiber
    
    Can contain the dispersion only (second derivative of the propagation constant) or a list of slopes from second order to nth order 
    
"""
struct Dispersion
    βₛ::Union{Float64,Vector{Float64}}
end

"""
    dispersion(D::Float64, λ::Float64)

    Generates the Dispersion structure for a fiber
    - D : dispersion factor in s/m²
    - λ : central wavelength of the study in m

"""
function dispersion(D::Float64, λ::Float64)
    Dispersion(-D * λ^2 / (2pi * c))
end