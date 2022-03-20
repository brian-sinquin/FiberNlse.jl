module FiberNlse
using FFTW
using ProgressBars

# Physical units & constants

nm = ns = 1e-9
ps = pm = 1e-12
km = 1e3
mW = mm = 1e-3
GHz = 1e9
Thz = 1e12
m = 1
W = 1
c = 299792458

struct NLSE_sim # Simulation struct containing all needed data
    Ψ::Matrix{Complex{Float32}}
    dz
    dt
    Nₜ::Int
    Nₗ::Int
    D
    γ
    α
    L
    T
    λ
    β2
    t::Vector{Float32}
    conf::SimulationConfig
end

struct Fiber
    D::Union{Float32} # dispersion factor
    α # exponential power loss
    γ # nonlinear Kerr factor
    L # fiber length
end

struct SimulationConfig
    SS::Bool
    SRS::Bool
    SBS::Bool
end

default_config = SimulationConfig(false, false, false)
function configure(Nₜ::Int,
    Nₗ::Int,
    fib::Fiber,
    T,
    λ, conf::SimulationConfig=default_config)

    return configure(Nₜ::Int,
        Nₗ::Int,
        fib.D,
        fib.γ,
        fib.α,
        fib.L,
        T,
        λ, conf)

end

function configure(
    Nₜ::Int,
    Nₗ::Int,
    D,
    γ,
    α,
    L,
    T,
    λ, conf::SimulationConfig=default_config)
    t = T * 0.5 * range(-1, stop=1, length=Nₜ) # Time vector
    Ψ = Matrix{Complex}(zeros((Nₗ, Nₜ)))
    return NLSE_sim(Ψ, L / Nₗ, T / Nₜ, Nₜ, Nₗ, D, γ, α, L, T, λ, -D * λ^2 / (2 * pi * c), t, conf), t
end

function inputSignal(sim, ψₒ)
    sim.Ψ[1, :] = ψₒ
end


function simulate(sim, progress::Bool=false)


    D(sim, ν) = @. -sim.β2 * 0.5im * (2 * pi * ν)^2
    N̂(sim, u) = @. abs(u)^2 * sim.γ * 1im

    # Check wether to show progressbar or not
    if progress
        iter_z = ProgressBar(range(2, sim.Nₗ))
    else
        iter_z = range(2, sim.Nₗ)
    end

    ν = FFTW.fftfreq(sim.Nₜ, 1.0 / (sim.dt))


    dl = sim.dz
    α = sim.α

    D̂ = D(sim, ν)

    for i in iter_z

        ψₗ = sim.Ψ[i-1, :]

        ψₗ = ifft(exp.(0.5 * dl .* (D̂ .- 0.5α)) .* fft(ψₗ))
        ψₗ = exp.(dl * N̂(sim, ψₗ)) .* ψₗ
        ψₗ = ifft(exp.(0.5 * dl .* (D̂ .- 0.5α)) .* fft(ψₗ))

        sim.Ψ[i, :] = ψₗ
    end

end

function transition(sim1, sim2) # Use the final state of sim1 as initial state of sim2
    inputSignal(sim2, sim1.Ψ[end, :])
end

end