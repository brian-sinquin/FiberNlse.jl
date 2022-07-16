"""
    Field(ψ::Matrix{ComplexF64}, l::Vector{Float64}, t::Vector{Float64})

    Structure containg space/time data resulting from the simulation
    - ψ : The signal propagation matrix 
    - l : The space vector in m 
    - t : The time vector in s
    
"""
struct Field
    ψ::Matrix{ComplexF64}
    l::Vector{Float64}
    t::Vector{Float64}
end




