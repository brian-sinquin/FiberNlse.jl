
<p align="center">
<img src="images/logo.png" />
</p>

---

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8251777.svg)](https://doi.org/10.5281/zenodo.8251777) [![CI](https://github.com/brian-sinquin/FiberNlse.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/brian-sinquin/FiberNlse.jl/actions/workflows/CI.yml) [![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/N/NamedDims.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html) [![codecov](https://codecov.io/gh/curio-sitas/FiberNlse.jl/branch/main/graph/badge.svg?token=O9L9P064J1)](https://codecov.io/gh/curio-sitas/FiberNlse.jl) [![code style blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle) [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

Documentation : [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://curio-sitas.github.io/FiberNlse.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://curio-sitas.github.io/FiberNlse.jl/dev)

---

 A non-linear Schrödinger equation solver for julia aimed towards fiber optics. 
**(🚧 Under development)**
## Features

The *FiberNlse.jl* package simulates the propagation of an optical field envelope signal of duration **T** in an optical fiber of length *L*. The chromatic dispersion (**D**) and SPM (self-phase modulation) (**γ**) wich arises from **Kerr** non-linearity are taken as parameters.

The core of the simulation consists in the integration of the *Non-Linear Schrödinger Equation* with the desired signal as initial condition. The package uses the Fourier Split-Step Method algorithm. 

---

## Installation
To install you can simply type :
```
] add FiberNlse
```
in your julia terminal or clone this repository and include the `src/FiberNlse.jl` file in your project.

## Roadmap

- [x] Implement Split-Step Method
- [x] Register DOI
- [x] Document code
- [x] Setup continuous integration
- [x] Add progress bar option
- [x] Add non constant dispersion (and higher order dispersion)
- [x] Add Self-steepening

     
- ### Future version (2.0) -- (compatibility breaking)
- [ ] Higher order integral solver (**DifferentialEquations.jl**)
- [ ] Add wavelength dependence to dispersion
- [ ] More complex material model API (GNLSE, Raman, SS, ...)
- [ ] Many solvers (as in [pychi](https://github.com/pychi-code/pychi/tree/main))
- [ ] Solver API (as **DifferentialEquations.jl**)
- [ ] Better Data manipulation (maybe a standalone package?)
- [ ] Make it compatible with **AD** for optimisation
      

## Citation

Please cite this repository if you use it to publish data in a research paper.


```
@software{sinquin_fibernlse,
  author       = {Sinquin Brian},
  title        = {FiberNlse.jl},
  month        = mar,
  year         = 2022,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.8251777},
  url          = {https://doi.org/10.5281/zenodo.8251777}
}
```

