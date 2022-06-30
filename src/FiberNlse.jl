module FiberNlse
export nm, ns, ps, pm, km, mW, mm, GHz, THz, m, W, c, smf28,edfa, simulate, configure, inputSignal, Fiber, Bundle
using FFTW
using ProgressBars


include("datatypes.jl")
# Physical units & constants

nm = ns = 1e-9
ps = pm = 1e-12
km = kHz = 1e3
mW = mm = 1e-3
GHz = 1e9
Thz = 1e12
m = 1
W = 1
c = 299792458





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
    l = L * range(0,stop=1,length=Nₗ)
    Ψ = Matrix{ComplexF32}(zeros((Nₗ, Nₜ)))
    return Simulation(Ψ, L / Nₗ, T / Nₜ, Nₜ, Nₗ, D, γ, α, L, T, λ, -D * λ^2 / (2 * pi * c), t, conf), t, l
end

function configure(
    Nₜ::Int,
    Nₗ::Int,
    fiber::Fiber,
    T,
    λ, conf::SimulationConfig=default_config)
    t = T * 0.5 * range(-1, stop=1, length=Nₜ) # Time vector
    l = fiber.L * range(0,stop=1,length=Nₗ)
    Ψ = Matrix{ComplexF32}(zeros((Nₗ, Nₜ)))
    return Simulation(Ψ, fiber.L / Nₗ, T / Nₜ, Nₜ, Nₗ, fiber.D, fiber.γ, fiber.α, fiber.L, T, λ, - fiber.D * λ^2 / (2 * pi * c), t, conf), t, l
end

function configure(
    fibs::Bundle,
    Nₜ::Int,
    T,
    λ, conf::SimulationConfig=default_config)


    t = T * 0.5 * range(-1, stop=1, length=Nₜ) # Time vector
    sims = Simulation[];

    l=Float64[];
    for (idx,fiber) in enumerate(fibs.elems)
        if typeof(fibs.Ns)==Int
            
                Nₗ = fibs.Ns;     
        else 
                
                Nₗ = fibs.Ns[idx];
        end
                
                li = fiber.L * range(0,stop=1,length=Nₗ)
                if isempty(l)
                    l=vcat(l,li)
                else
                    l=vcat(l,li.+l[end])
                end
                Ψ = Matrix{ComplexF32}(zeros((Nₗ, Nₜ)))
                sim = Simulation(Ψ, fiber.L / Nₗ, T / Nₜ, Nₜ, Nₗ, fiber.D, fiber.γ, fiber.α, fiber.L, T, λ, - fiber.D * λ^2 / (2 * pi * c), t, conf)
        push!(sims,sim)
    end


   
    return sims,t,l
end
function inputSignal(sim, ψₒ)
    sim.Ψ[1, :] = ψₒ
end

function edfa(L::Real, G::Real)
    g  = 0.1*log(10)*G/L
    fib = Fiber(-10*ps/nm/km, -g, 1/W/km, L);
    return fib
end

function smf28(L::Real)
    fib = Fiber(17*ps/nm/km, 0.046/km, 1.1/W/km, L);
    return fib
end



function simulate(sim::Simulation, progress::Bool=false)


    D(sim, ν) = @. sim.β2 * 0.5im * (2 * pi * ν)^2
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

function simulate(sims::Vector{Simulation}, Ψₒ::Vector{Float64},progress::Bool=false)
    S = []
    for (idx, sim) in enumerate(sims)
        if idx == 1
            inputSignal(sim,Ψₒ)
            simulate(sim,progress)
            S = sim.Ψ
        else
            inputSignal(sim, S[end, :])
            simulate(sim,progress)
            S = vcat(S, sim.Ψ)
        end
        
    end
    return S
end

function transition(sim1::Simulation, sim2::Simulation) # Use the final state of sim1 as initial state of sim2
    inputSignal(sim2, sim1.Ψ[end, :])
end


function rhs()
    # ! Lire le code de gnlse-python (les mettre dans les crédits)
end

function simulate2(sim::Simulation)
    # ! Integration RK4 de l'ode rhs
end
end