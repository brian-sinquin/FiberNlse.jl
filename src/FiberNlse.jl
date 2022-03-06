module FiberNlse
using FFTW
using ProgressBars


struct NLSE_sim # Simulation struct containing all needed data
   Ψ::Matrix{Complex{Float64}}
   dz::Float32
   dt::Float32
   N_t::Int
   N_z::Int
   D::Float64
   γ::Float64
   α::Float64
   L::Int
   T::Float64
   λ::Float64
   β2::Float64
end

function configure(
    N_t::Int,
    N_z::Int,
    D::Float64,
    γ::Float64,
    α::Float64,
    L::Float64,
    T::Float64,
    λ::Float64)

    Ψ = Matrix{Complex{Float64}}(zeros((N_z, N_t)))
    return NLSE_sim(Ψ, L/N_z, T/N_t, N_t, N_z, D, γ, α, L, T, λ,-D*λ^2/(2*pi*c))
end

function initialSignal(sim, ψₒ::Vector{ComplexF32})
    sim.Ψ[1,:] = ψₒ;
end

function simulate(sim, progress::Bool=false)

    
    ν = FFTW.fftfreq(sim.N_t, 1.0/(2*sim.dt))

    Disp(β₂) = @. -β₂*0.5im*(2*pi*ν)^2
    N(ui) = abs.(ui).^2*sim.γ*1im
    dl = sim.dz
    α = sim.α
    if progress
        iter_z = ProgressBar(range(2, sim.N_z));
    else 
        iter_z = range(2, sim.N_z);
    end
    for i in iter_z

        ψₗ=sim.Ψ[i-1,:];

        ψₗ = ifft(exp.(0.5*dl.*(Disp(sim.β2).-0.5.*α)).*fft(ψₗ));
        ψₗ = exp.(dl*N(ψₗ)).*ψₗ;
        ψₗ = ifft(exp.(0.5*dl.*(Disp(sim.β2).-0.5.*α)).*fft(ψₗ));

        sim.Ψ[i,:] = ψₗ;
    end

end

function branch(sim1, sim2) # Use the final state of sim1 as initial state of sim2
    sim2.initialSignal(sim1.Ψ[end,:]);
end

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



export Units

end