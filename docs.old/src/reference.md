# Technical Reference

This section provides detailed technical information about parameters and characteristic scales used in nonlinear fiber optics and waveguide simulations.

## Unit Conventions

FiberNlse.jl uses SI units throughout, following conventions established in standard nonlinear fiber optics literature[^Agrawal][^Dudley].

### Base Units and Typical Ranges

| Parameter | Symbol | SI Unit | Typical Values | Notes |
|-----------|--------|---------|----------------|-------|
| Time | t | s | 10⁻¹⁵ - 10⁻¹² | Femtosecond to picosecond pulses |
| Distance | z | m | 10⁻³ - 10³ | mm to km propagation distances |
| Power | P | W | 1 - 10⁴ | W to kW peak powers |
| Wavelength | λ | m | ~1.55×10⁻⁶ | Typically C-band for telecom |

### Dispersion Parameters[^Agrawal_ch3]

| Parameter | Symbol | SI Unit | Glass Fiber | Silicon | ChG |
|-----------|--------|---------|-------------|---------|-----|
| GVD | β₂ | s²/m | ±20×10⁻²⁷ | -0.9×10⁻²⁴ | -400×10⁻²⁷ |
| TOD | β₃ | s³/m | ~10⁻⁴⁰ | ~10⁻³⁹ | ~10⁻⁴⁰ |
| D | D | s/m² | ~17×10⁻⁶ | ~700×10⁻⁶ | ~300×10⁻⁶ |

The dispersion parameter D (commonly given in ps/nm/km) relates to β₂ as:

```math
D = -\frac{2\pi c}{\lambda^2}\beta_2
```

### Nonlinear Parameters[^Agrawal_ch2]

| Parameter | Symbol | SI Unit | Glass Fiber | Silicon | ChG |
|-----------|--------|---------|-------------|---------|-----|
| Nonlinearity | γ | W⁻¹m⁻¹ | ~10⁻³ | 100-500 | 10-20 |
| Loss | α | m⁻¹ | ~10⁻⁴ | 100-200 | ~1 |
| Raman fraction | fᵣ | - | 0.18 | - | 0.03 |

## Characteristic Lengths[^Agrawal_ch3]

For proper numerical modeling, it's important to consider these characteristic length scales:

### Dispersion Length
```math
L_D = \frac{T_0^2}{|\beta_2|}
```
This length characterizes the distance over which dispersive effects become important. When ``z \sim L_D``, the pulse shape begins to change significantly due to dispersion.

### Nonlinear Length
```math
L_{NL} = \frac{1}{\gamma P_0}
```
This length characterizes the distance over which nonlinear effects become important. When ``z \sim L_{NL}``, the pulse spectrum begins to broaden due to self-phase modulation.

### Walk-off Length
```math
L_W = \frac{T_0}{|\beta_1(\omega_1) - \beta_1(\omega_2)|}
```
This length is relevant for multi-wavelength interactions, characterizing the distance over which pulses at different wavelengths separate due to group velocity mismatch.

## Numerical Considerations

The characteristic lengths help determine appropriate simulation parameters:

1. Step size selection: Should be much smaller than all relevant characteristic lengths
2. Total propagation distance: Compare to the shortest characteristic length
3. Grid resolution: Should resolve the shortest temporal features
4. Window size: Should accommodate temporal walk-off and spectral broadening

[^Agrawal]: G. P. Agrawal, "Nonlinear Fiber Optics," Academic Press, 6th Edition (2019).
[^Agrawal_ch2]: Chapter 2 of Agrawal's "Nonlinear Fiber Optics" covers nonlinear effects and their parameters in detail.
[^Agrawal_ch3]: Chapter 3 of Agrawal's "Nonlinear Fiber Optics" provides comprehensive coverage of dispersion and characteristic lengths.
[^Dudley]: J. M. Dudley et al., "Supercontinuum generation in photonic crystal fiber," Reviews of Modern Physics 78, 1135 (2006).
