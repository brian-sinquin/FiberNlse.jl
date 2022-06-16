struct Fiber
    D::Union{Float32} # dispersion factor
    α # exponential power loss
    γ # nonlinear Kerr factor
    L # fiber length
end

struct Bundle
    elems::Vector{Fiber}
    Ns::Union{Int, Vector{Int}} 
end

struct SimulationConfig
    SS::Bool
    SRS::Bool
    SBS::Bool
end

default_config = SimulationConfig(false, false, false)

struct Simulation # Simulation struct containing all needed data
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
