function propagate(ψ₀::Vector{ComplexF64}, fib::Fiber, T::Float64, Nₗ::Int)
    Nₜ=length(ψ₀)
    dt,dz = T/Nₜ, fib.L/Nₗ

    t = (0:Nₜ-1)dt
    l = (0:Nₗ-1)dz

    α = fib.α
    γ = fib.γ
    ν = FFTW.fftfreq(Nₜ, 1. /dt)

    if typeof(fib.D.βₛ)==Float64 D̂ = @. fib.D.βₛ * 2im * ( pi * ν)^2 end

    N̂(u) = @. abs(u)^2 * γ * im

    # Check wether to show progressbar or not
    #if progress
        #iter_z = ProgressBar(range(2, sim.Nₗ))
    #else
        #iter_z = range(2, sim.Nₗ)
    #end

    ψ = Matrix{ComplexF64}(zeros((Nₗ, Nₜ)))
    ψ[1,:] = ψ₀

    for i in 2:Nₗ

        ψ[i, :] = ifft(exp.(0.5 * dz .* (D̂ .- 0.5α)) .* fft(ψ[i-1, :]))
        ψ[i, :] = exp.(dz * N̂(ψ[i, :])) .* ψ[i, :]
        ψ[i, :] = ifft(exp.(0.5 * dz .* (D̂ .- 0.5α)) .* fft(ψ[i, :]))

    end
    Field(ψ, l, t)

end