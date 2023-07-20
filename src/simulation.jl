using ProgressBars

# Future SSFM(...)
function propagate(
    ψ₀::Union{Vector{ComplexF64},Vector{Float64}},
    fib::Fiber,
    T::Float64,
    Nₗ::Int;
    progress = false,
    convention::Symbol = getPhaseConvention(),
)
    Λ = getPhaseConventionSign(convention)
    Nₜ = length(ψ₀)
    dt, dz = T / Nₜ, fib.L / Nₗ

    t = (-Nₜ÷2:Nₜ÷2-1)dt
    l = (0:(Nₗ-1))dz

    α = fib.α
    γ = fib.γ
    ν = fftshift(collect(-Nₜ/2:Nₜ/2-1)/(Nₜ*dt)); #frequency grid
    ω = 2pi*ν; #frequency grid

    # Dispersion vector from propagation constants Taylor expansion around λ
    D̂ = Λ .* sum(fib.D.β[i] .* (im)^(i) .* (1im .* ω) .^ (i + 1) ./ factorial(i + 1) for i in 1:length(fib.D.β))

    # Nonlinear operator including self-steepening
    N̂(u) = - Λ .* γ * im .* (abs.(u) .^ 2 .- ifft(ν .* fft(u .* abs.(u) .^ 2)) .* fib.λ / c)

    # Check wether to show progressbar or not
    if progress
        iter = ProgressBar(2:Nₗ)
    else
        iter = 2:Nₗ
    end

    ψ = zeros(ComplexF64,Nₗ, Nₜ)
    ψ[1, :] = complex(ψ₀)

    for i in iter
        #TODO enhance matrix allocation (plan_fft)
        ψ[i, :] = ifft(exp.(0.5 * dz .* (D̂ .- 0.5α)) .* fft(@view ψ[i-1, :]))
        ψ[i, :] = exp.(dz * N̂(@view ψ[i, :])) .* (@view ψ[i, :])
        ψ[i, :] = ifft(exp.(0.5 * dz .* (D̂ .- 0.5α)) .* fft(@view ψ[i, :]))
    end
    return Field(ψ, l, t)
end

# Future SSFM(...)
function propagate(
    ψ₀::Union{Vector{ComplexF64},Vector{Float64}},
    fibs::Vector{Fiber},
    T::Float64,
    Nₗ::Int;
    progress = false,
    convention::Symbol = getPhaseConvention(),
)
    field = propagate(ψ₀, fibs[1], T, Nₗ; progress, convention)
    for i = 2:length(fibs)
        #TODO enhance matrix allocation
        field = concatf(field, propagate(output(field), fibs[i], T, Nₗ; progress))
    end
    return field
end

# Add GNLSE 
# Find a common structure for time/freq field object (common for both methods)