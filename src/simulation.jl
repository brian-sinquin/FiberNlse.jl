using ProgressBars

function propagate(ψ₀::Union{Vector{ComplexF64},Vector{Float64}}, fib::Fiber, T::Float64, Nₗ::Int; progress=false)
    Nₜ=length(ψ₀)
    dt,dz = T/Nₜ, fib.L/Nₗ

    t = (0:Nₜ-1)dt
    l = (0:Nₗ-1)dz

    α = fib.α
    γ = fib.γ
    ν = FFTW.fftfreq(Nₜ, 1. /dt)

    # Dispersion vector from propagation constants Taylor expansion around λ
    if typeof(fib.D.β)==Float64
        D̂ = @. fib.D.β * 2im * ( pi * ν)^2
    else
        D̂ = zeros(length(ν)) .+ 0 .*im
        for i in 1:length(fib.D.β)
            D̂ = D̂ .- fib.D.β[i] .* (im)^(i) .* (2pi.*im.* ν).^(i+1)./factorial(i+1)
        end
    end


    # Nonlinear operator including self-steepening
    N̂(u) =  γ*im .* (abs.(u).^2 .-ifft(ν.*fft(u.*abs.(u).^2)).*fib.λ/c)

    # Check wether to show progressbar or not
    if progress
        iter = ProgressBar(2:Nₗ)
    else
        iter = 2:Nₗ
    end

    ψ = Matrix{ComplexF64}(zeros((Nₗ, Nₜ)))
    ψ[1,:] = complex(ψ₀)

    for i in iter
        #TODO enhance matrix allocation (plan_fft)
        ψ[i, :] = ifft(exp.(0.5 * dz .* (D̂ .- 0.5α)) .* fft(@view ψ[i-1, :]))
        ψ[i, :] = exp.(dz * N̂(@view ψ[i, :])) .* (@view ψ[i, :])
        ψ[i, :] = ifft(exp.(0.5 * dz .* (D̂ .- 0.5α)) .* fft(@view ψ[i, :]))

    end
    Field(ψ, l, t)

end

function propagate(ψ₀::Union{Vector{ComplexF64},Vector{Float64}}, fibs::Vector{Fiber}, T::Float64, Nₗ::Int; progress=false)
   field = propagate(ψ₀, fibs[1], T, Nₗ; progress)
   for i in 2:length(fibs)
        #TODO enhance matrix allocation
        field = concatf(field,propagate(output(field), fibs[i], T, Nₗ; progress))
   end
   field
end
