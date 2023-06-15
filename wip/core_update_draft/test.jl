using Revise
using Plots

include("gnlse.jl")

#numerical grid
n = 2^13;                   #number of grid points
twidth = 200e-12;          #width of time window [s]
c = 299792458;              #speed of light [m/s]
wavelength = 1550e-9;        #reference wavelength [m]

dt = twidth/n;
T = (-n/2:n/2 - 1).*dt; #time grid

# input pulse
power = 1;              #peak power of input [W]
t0 = 30e-12;              #duration of input [s]
C=-0.0
A = @. sqrt(power)*exp(-(T/t0)^2*(1+1im*C)); #input field [W^(1/2)]

# fibre parameters
flength = 5000.0;             #fibre length [m]
#betas = [beta2, beta3, ...] in units [s^2/m, s^3/m ...]
betas = [-1.0e-026,0.3e-038]
gamma = 1.1e-3;               #nonlinear coefficient [1/W/m]
loss = 0.2e-3;                   #loss [dB/m]

#simulation parameters
nsaves = 200;     #number of length steps to save field at


#propagate field
Z, AT, AW, W = gnlse(T, A, w0, gamma, betas, loss, flength, nsaves; raman=true);
                   
# plot output

lims = (1540, 1570)
WL = reverse(2pi*c./W)
iis = findall(lims[1] .<= WL*1e9 .<= lims[2])
dBIW = reverse(10log10.(abs2.(AW) .* 2pi*c./WL'.^2), dims=2)
dBIWmax = maximum(dBIW)
heatmap(WL[iis]*1e9, Z, dBIW[:,iis], clim=(dBIWmax-40.0, dBIWmax))


dBIT = 10log10.(abs2.(AT))
dBITmax = maximum(dBIT)

heatmap(T*1e12, Z, abs2.(AT))

plot(abs2.(AT[end,:]))
plot!(abs2.(AT[1,:]))

function derivate(y::AbstractVector, x::AbstractVector)
    function centraldiff(v::AbstractVector)
        dv = diff(v) / 2 # half the derivative
        a = [dv[1]; dv] # copies first element
        a .+= [dv; dv[end]] # copies last element, add both results to compute average
        return (a)
    end
    return centraldiff(y) ./ centraldiff(x)
end
using DSP

phase = unwrap(angle.(AT), dims=2)
plot(phase[end,:])
heatmap(phase, cm=:spectrum)


@gif for i ∈ 1:nsaves
    plot(T*1e12,abs2.(AT[i,:]), ylims=(0,20), color=:black)
    plot!(twinx(),T*1e12,unwrap(angle.(AT[i,:])),ylims=(0,20), color=:red, ls=:dash, label="z=$(round(Int,Z[i])) m")
    xlims!(-50,50)
end 