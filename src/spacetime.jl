"""
    Field(ψ::Matrix{ComplexF64}, l::Vector{Float64}, t::Vector{Float64})

    Structure containing space/time data resulting from the simulation
    - ψ : The signal propagation matrix
    - l : The space vector in m
    - t : The time vector in s

"""
struct Field
    ψ::Matrix{ComplexF64}
    l::Vector{Float64}
    t::Vector{Float64}
end

"""
    output(f::Field)

    returns the signal at the end of the propagation
"""
function output(f::Field) f.ψ[end,:] end

#(f1::Field,f2::Field) = concatf(f1,f2)

function concatf(f1::Field, f2::Field)
    return Field(vcat(f1.ψ, f2.ψ), vcat(f1.l, f2.l .+ f1.l[end]), f2.t)
end
