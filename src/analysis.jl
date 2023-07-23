
function phase(ψ; convention::Symbol = getPhaseConvention())
    sign = getPhaseConventionSign(convention)
    sign * (ψ .|> angle |> unwrap)
end

function instFreq(ψ, t, convention::Symbol = getPhaseConvention())
    derivate(phase(ψ; convention), t) / (2pi)
end

function spectrum(ψ::Vector{ComplexF64}, dt)
    N = length(ψ);
    ν = collect(-N/2:N/2-1)/(N*dt);
    sp = fftshift(fft(ψ)/N)
    (ν,sp)
end
