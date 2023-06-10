using FiberNlse, Plots 

c=3e8

λ = 1.55e-6



begin # smf
Lsmf = 1000.0 
Dsmf = 18ps/nm/km
Ssmf = 0.09ps/(km*nm^2)
smf_β₂ = -λ^2*Dsmf/(2pi*c)
smf_β₃ = -(λ^2/(2pi*c))^2*Ssmf
smf = Fiber(Lsmf, Dispersion([smf_β₂,smf_β₃]), 1.1e-3, log(10)*0.2/1000/10, λ)
end

begin # hnlf
    Lhnlf = 100.0 
    Dhnlf = 0.46ps/nm/km
    Shnlf = 0.018ps/(km*nm^2)
    hnlf_β₂ = -λ^2*Dhnlf/(2pi*c)
    hnlf_β₃ = -(λ^2/(2pi*c))^2*Shnlf
    hnlf = Fiber(Lhnlf, Dispersion([hnlf_β₂,hnlf_β₃]), 10.1e-3, log(10)*0.8/1000/10, λ)
end

begin
n = 2^15
T = 200e-12
t = 0.5*T*LinRange(-1,1, n)
τ = 25ps
Pc0 = 6.1
E0 = sqrt(Pc0)*exp.(-(t/τ).^2)

field1 = propagate(E0, smf, T, 1000)
Eout1 = output(field1)

Pc1 = 7#maximum(abs2.(Eout1))
field2 = propagate(sqrt(Pc1)*Eout1/sqrt(maximum(abs2.(Eout1))), hnlf, T, 1000)
Eout2 = output(field2)
end

begin
plot(t/ps, E0 .|> abs2)
plot!(t/ps, Eout1 .|> abs2)
plot!(t/ps, Eout2 .|> abs2)
xlims!(-100,100)
end

idt = findall(abs.(t).<= 10ps)

M = field2.ψ[:,idt] .|> abs2
heatmap(M)
using DSP: fft, fftshift
function spectrum(s,t)
    dt = t[2]-t[1]
    f = LinRange(-1,1,length(t))/(2*dt)
    sp = fftshift(fft(s))/length(s)
    return sp,f
end

normalize(x)=x./maximum(abs.(x)) 
logsp(x) = 10log10.(x)

begin
sp0, f = spectrum(Eout1,t)
SP0 = sp0 .|> abs2  |> normalize  |> logsp
sp1, _ = spectrum(Eout2,t)
SP1 = sp1 .|> abs2  |> normalize  |> logsp

plot(f/1e12, SP0)
plot!(f/1e12, SP1)
xlims!(-20,20)
ylims!(-40,0)
end