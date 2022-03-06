<p align="center">
<img src="logo.png" />
</p>

# 
 A non-linear Schrödinger equation solver for julia aimed towards fiber optics.

## Features

The *FiberNlse.jl* package simulates the propagation of an optical field envelope signal of duration *T* in an optical fiber of length *L*. The chromatic dispersion *D* and SPM (self-phase modulation) *γ* wich arises from **Kerr** non-linearity are taken as parameters.

The core of the simulation consists in the integration of the *Non-Linear Schrödinger Equation* with the desired signal as inital condition.

## Examples
- [ ] TODO
### Sech² Soliton
### Self-phase modulation
### Time-lens compression

## Installation
To install you can simply type `] add FiberNlse` in your julia terminal or clone this repository and include the `src/FiberNlse.jl` file in your project.

## Roadmap

- [x] Implement Split-Step Method
- [x] Add progress bar option
- [ ] Document code
- [ ] Setup continuous integration
- [ ] Add non constant dispersion (and higher order dispersion)
- [ ] Create a documentation page (mdbook ?)
- [ ] Higher order integral solver (**DifferentialEquations.jl**)
- [ ] Add more non-linear processes (Raman & Brillouin scattering)

### source : https://www.fiberoptics4sale.com/blogs/wave-optics/solitons-in-optical-fibers