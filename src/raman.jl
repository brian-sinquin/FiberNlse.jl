"""Raman scattering function for silica optical fibers, based on Q. Lin
and Govind P. Agrawal model.
Parameters
 ----------
T : float
    Time vector.
Returns
-------
fr : float
   Share of Raman response.
RT : ndarray
   Vector representing Raman response.
"""
function raman_linagrawal(T)

    # Raman response [arbitrary units]
    fr = 0.245
    # Adjustable parameters used to fit the actual Raman gain spectrum [ps]
    tau1 = 0.0122
    tau2 = 0.032
    taub = 0.096
    # Fractional contribution of the anisotropic reponse to the total Raman
    # response
    fb = 0.21
    fc = 0.04
    # Fractional contribution of the isotropic reponse to the total Raman
    # response
    fa = 1 - fb - fc
    # Anisotropic Raman response
    ha = @. (tau1^2 + tau2^2) / tau1 / (tau2^2) * exp(-T / tau2) * sin(T / tau1)
    # Izotropic Raman respons
    hb = @. (2 * taub - T) / (taub^2) * exp(-T / taub)
    # Total Raman response
    RT = @. (fa + fc) * ha + fb * hb

    RT[findall(T .< 0)] .= 0

    return fr, RT
end
