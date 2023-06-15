using FFTW
using Trapz
using DifferentialEquations
using ProgressMeter


## begin function
function gnlse(T, A, w0, gamma, betas, loss, flength, nsaves; progress=true, raman=false)

  if raman 
  fr = 0.18;                  #fractional Raman contribution
  tau1 = 0.0122e-12; tau2 = 0.032e-12;
  RT = @. (tau1^2+tau2^2)/tau1/tau2^2*exp(-T/tau2).*sin(T/tau1);
  RT[findall(T .< 0)] .= 0;                #heaviside step function
  RT = RT/(trapz(T, RT))
  else
    fr = 0
    RT = zeros(length(T))
  end

n = length(T); dT = T[2]-T[1]; #grid parameters
V = 2*pi*collect(-n/2:n/2-1)/(n*dT); #frequency grid
alpha = log(10 .^(loss/10));    #attenuation coefficient

B = zeros(Float64, length(V));
for i = 1:length(betas)        #Taylor expansion of betas
B = B .+ betas[i]/factorial(i+1).*V.^(i+1);
end
L = 1im*B .- alpha/2;            #linear operator

if abs(w0) > eps()             #if w0>0 then include shock
gamma = gamma/w0;    
W = V .+ w0;                #for shock W is true freq
else
W = ones(length(V));                     #set W to 1 when no shock
end

RW = n*ifft(fftshift(RT));   #frequency domain Raman
L = fftshift(L); W = fftshift(W); #shift to fft space

#define function to return the RHS of Eq. (3.13)
p = Progress(100)
function progress_callback(int)
    update!(p, round(Int, int.t/flength*100))
end

cb = PeriodicCallback(progress_callback, flength/100)
function rhs!(dAW, AW, p, z) 

AT = fft(AW.*exp.(L.*z));         #time domain field
IT = abs.(AT).^2;                #time domain intensity
if (length(RT) == 1) || (abs(fr) < eps()) #no Raman case
M = ifft(AT.*IT);             #response function
else
RS = dT*fr*fft(ifft(IT).*RW); #Raman convolution
M = ifft(AT.*((1-fr).*IT .+ RS)); #response function
end
dAW[:] = 1im*gamma*W.*M.*exp.(-L.*z);   #full RHS of Eq. (3.13)
end


# setup and run the ODE integrator
Z = collect(LinRange(0.0, flength, nsaves));  #select output z points
# set error control options
#options = odeset('RelTol', 1e-5, 'AbsTol', 1e-12, ...
#'NormControl', 'on', ...
#'OutputFcn', @report);
problem = ODEProblem(rhs!, ifft(A), (Z[1],Z[end]));
sol = solve(problem, alg=DP5(), reltol=1e-6, abstol=1e-6, saveat=(Z[2]-Z[1]), callback = cb)
#sol = ode45(rhs, Z, ifft(A), options); #run integrator
AW = zeros(ComplexF64, nsaves, length(A))
for i in 1:nsaves
    AW[i,:] = sol(Z[i]);
end
# process output of integrator
AT = zeros(ComplexF64, size(AW));
for i = 1:size(AW,1)
  AW[i,:] = AW[i,:].*exp.(L.*Z[i]); # change variables
  AT[i,:] = fft(AW[i,:]);           # time domain output
  AW[i,:] = fftshift(AW[i,:]).*dT.*n;  # scale
end

W = V .+ w0; #the absolute frequency grid

return (Z, AT, AW, W)
end