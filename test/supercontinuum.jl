using FiberNlse, Plots, StatsBase, DSP
FiberNlse.setPhaseConvention(:positive)
FiberNlse.getPhaseConvention()

#numerical grid
n = 2^13;                   #number of grid points
twidth = 12.5e-12;          #width of time window [s]
c = 299792458;              #speed of light [m/s]
wavelength = 835e-9;        #reference wavelength [m]

dt = twidth/n;
T = (-n/2:n/2 - 1).*dt; #time grid

# input pulse
power = 10000;              #peak power of input [W]
t0 = 28.4e-15;              #duration of input [s]
A =  @. sqrt(power)*sech(T/t0); #input field [W^(1/2)]

# fibre parameters
flength = 0.15;             #fibre length [m]
#betas = [beta2, beta3, ...] in units [s^2/m, s^3/m ...]
betas = [-1.1830e-026, 8.1038e-041, -9.5205e-056,  2.0737e-070, 
         -5.3943e-085,  1.3486e-099, -2.5495e-114,  3.0524e-129,
         -1.7140e-144];
gamma = 0.11;               #nonlinear coefficient [1/W/m]
loss = 0.0;                   #loss [dB/m]

#simulation parameters
nsaves = 200;     #number of length steps to save field at

fib = Fiber(flength, Dispersion(betas),  gamma, loss, wavelength)
#propagate field

field = propagate(A, fib, twidth, nsaves; progress=true, raman=true)

plot(10log10.((field.ψ[1,:]) .|> abs))
plot!(10log10.((field.ψ[end,:]) .|> abs))
using FFTW
AW = field.AW

WL =reverse(2pi*c./(field.W))

lims = (450, 1350)

iis = findall(lims[1] .<= WL*1e9 .<= lims[2])
dBIW = reverse(10log10.(abs2.(AW) .* 2pi*c./WL'.^2), dims=2)
dBIWmax = maximum(dBIW)
heatmap(WL[iis]*1e9, field.l, dBIW[:,iis], clim=(dBIWmax-40.0, dBIWmax))


dbIt = 10log10.(abs2.(field.ψ))
heatmap(field.t*1e12, field.l, dbIt,dims=2, clim=(maximum(dbIt)-40,maximum(dbIt)))
xlims!(5.75,11.25)

ATend = field.ψ[end,:]
Awtend = fftshift(fft(ATend))
Awend = AW[end,:]

#! TODO scale spectrum ?? 
plot(WL*1e9,10log10.(abs2.(Awtend)))
plot!(WL*1e9,10log10.(abs2.(Awend) .* 2pi*c./WL.^2))
xlims!(450,3350)