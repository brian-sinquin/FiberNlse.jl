# FiberNlse.jl

## What is it ?

It is a **Fourier Split-Step** based algorithm which solves the **Non-Linear Schrödinger Equation**(NLSE).
It uses as parameters the usual quantities describing optical fibers and optical signal.
## To do what ?

The package is designed for **photonic researchers** using optical fibers who need to simulate the propagation of an optical signal **in such a medium**. It can be used to study the evolution of a **pulse train** or the degradation of a **transmission line** by dispersion and nonlinearities.
## Features

The package currently supports :
    * Kerr nonlinearity
    * Linear losses/gain
    * Optical envelope (SVEA)
    * Self-steepening
    * Dispersion (from order 2 to n)
## How to use it ?

You will learn to use the packages by reading the **user guide** and the **examples**
## Roadmap 

In the future, the plan is to implement :
    * Raman scattering (intra-pulse)
    * Brillouin scattering
    * Polarisation (PMD)
    * Random sequence generator and BER tests/ Eye diagram
