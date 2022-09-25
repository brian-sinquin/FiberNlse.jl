
function phase(ψ; convention::Symbol = getPhaseConvention())
    sign = getPhaseConventionSign(convention)
    sign * (ψ .|> angle |> unwrap)
end

function instFreq(ψ, t, convention::Symbol = getPhaseConvention())
    derivate(phase(ψ; convention), t) / (2pi)
end

function spectrum(ψ::Vector{ComplexF64}, dt)
    return error("To implement")
end
