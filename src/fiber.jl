"""
    Fiber(L::Float64, D::Dispersion, γ::Float64, α::Float64)

    Structure containing fiber caracteristics :
    - L : length in m 
    - D : dispersion
    - γ : nonlinear factor in 1 / W / m
    - α : fiber attenuation or gain (exponential) in 1 / m
"""
struct Fiber 
    L::Float64
    D::Dispersion
    γ::Float64
    α::Float64
end

function edfa(L::Float64, G::Float64, λ::Float64)
    g  = 0.1*log(10)*G/L
    fib = Fiber(L,disperion(-10e-6,λ), 1.1e-3, -g);
    return fib
end

function smf28(L::Float64, λ::Float64)
    fib = Fiber(L,dispersion(17e-6, λ), 1.1e-3, 0.046e-3);
    return fib
end