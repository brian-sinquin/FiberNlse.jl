function propagate(ψ₀::Union{Vector{ComplexF64},Vector{Float64}}, fib::Fiber, T::Float64, Nₗ::Int)
    Nₜ=length(ψ₀)
    dt,dz = T/Nₜ, fib.L/Nₗ

    t = (0:Nₜ-1)dt
    l = (0:Nₗ-1)dz

    α = fib.α
    γ = fib.γ
    ν = FFTW.fftfreq(Nₜ, 1. /dt)

    if typeof(fib.D.β)==Float64 D̂ = @. fib.D.β * 2im * ( pi * ν)^2 end

    N̂(u) = @. abs(u)^2 * γ * im

    # Check wether to show progressbar or not
    #if progress
        #iter_z = ProgressBar(range(2, sim.Nₗ))
    #else
        #iter_z = range(2, sim.Nₗ)
    #end

    ψ = Matrix{ComplexF64}(zeros((Nₗ, Nₜ)))
    ψ[1,:] = complex(ψ₀)
    #= self.fr, RT = setup.raman_model(self.t)
    self.RW = self.N * np.fft.ifft(np.fft.fftshift(np.transpose(RT)))
    IT = np.abs(At)**2

    if self.RW is not None:
        X[:] = IT
        plan_inverse()
        x[:] *= self.RW
        plan_forward()
        RS = dt * self.fr * X
        X[:] = At * ((1 - self.fr) * IT + RS)
        M = plan_inverse()
    else:
        X[:] = At * IT
        M = plan_inverse() =#

    for i in 2:Nₗ

        ψ[i, :] = ifft(exp.(0.5 * dz .* (D̂ .- 0.5α)) .* fft(ψ[i-1, :]))
        ψ[i, :] = exp.(dz * N̂(ψ[i, :])) .* ψ[i, :]
        ψ[i, :] = ifft(exp.(0.5 * dz .* (D̂ .- 0.5α)) .* fft(ψ[i, :]))

    end
    Field(ψ, l, t)

end

function propagate(ψ₀::Union{Vector{ComplexF64},Vector{Float64}}, fibs::Vector{Fiber}, T::Float64, Nₗ::Int)
   field = propagate(ψ₀, fibs[1], T, Nₗ)
   for i in 2:length(fibs)
        field = concatf(field,propagate(output(field), fibs[i], T, Nₗ))
   end
   field
end
