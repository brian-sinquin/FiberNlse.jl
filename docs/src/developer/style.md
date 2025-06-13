# Code Style Guide

FiberNlse.jl follows the [BlueStyle](https://github.com/invenia/BlueStyle) convention for Julia code. This guide outlines the key style points and best practices.

## General Guidelines

### Naming Conventions
- Use clear, descriptive names
- CamelCase for types and modules
- snake_case for functions and variables
- ALL_CAPS for constants

```julia
# Good
struct WaveguideParameters
    group_velocity_dispersion::Float64
    nonlinear_coefficient::Float64
end

# Bad
struct waveguide_params
    gvd::Float64
    gamma::Float64
end
```

### Function Naming
- Use verb-noun format for functions that perform actions
- Use noun format for functions that return values

```julia
# Good
calculate_dispersion(β::Vector{Float64}, ω::Vector{Float64})
propagate_field(field::Vector{ComplexF64}, distance::Float64)

# Bad
disp(β::Vector{Float64}, ω::Vector{Float64})
field_prop(field::Vector{ComplexF64}, distance::Float64)
```

## Code Organization

### File Structure
```julia
# 1. Module definition
module MyModule

# 2. Imports
using Package1
using Package2: specific_function

# 3. Exports
export public_function, PublicType

# 4. Constants
const SPEED_OF_LIGHT = 299792458.0

# 5. Types
struct MyType
    field1::Float64
    field2::Int
end

# 6. Functions
function my_function(x::Float64)
    # Implementation
end

end # module
```

### Type Hierarchy
Organize types from most abstract to most concrete:

```julia
abstract type AbstractWaveguide end

struct OpticalFiber <: AbstractWaveguide
    # Implementation
end

struct PhotonicCrystalFiber <: OpticalFiber
    # Implementation
end
```

## Documentation

### Function Documentation
Use complete docstrings:

```julia
"""
    calculate_dispersion(β::Vector{Float64}, ω::Vector{Float64}) -> Vector{Float64}

Calculate the dispersion operator in the frequency domain.

# Arguments
- `β::Vector{Float64}`: Dispersion coefficients
- `ω::Vector{Float64}`: Angular frequency grid

# Returns
- `Vector{Float64}`: Frequency-domain dispersion operator

# Example
```julia
β = [0.0, -20e-27, 1e-40]
ω = range(-1e15, 1e15, length=1024)
D = calculate_dispersion(β, ω)
```
"""
function calculate_dispersion(β::Vector{Float64}, ω::Vector{Float64})
    # Implementation
end
```

### Type Documentation
Document types with their purpose and fields:

```julia
"""
    Waveguide(loss, betas, gamma, wavelength, length)

A type representing an optical waveguide with its physical parameters.

# Fields
- `loss::Float64`: Linear power loss coefficient [dB/m]
- `betas::Vector{Float64}`: Dispersion coefficients [s^n/m]
- `gamma::Float64`: Nonlinear coefficient [1/W/m]
- `wavelength::Float64`: Reference wavelength [m]
- `length::Float64`: Waveguide length [m]
"""
struct Waveguide
    loss::Float64
    betas::Vector{Float64}
    gamma::Float64
    wavelength::Float64
    length::Float64
end
```

## Testing

### Test Organization
Organize tests by functionality:

```julia
# test/runtests.jl
using Test
using FiberNlse

@testset "FiberNlse.jl" begin
    @testset "Dispersion" begin
        # Dispersion-related tests
    end
    
    @testset "Nonlinearity" begin
        # Nonlinearity-related tests
    end
    
    @testset "Integration" begin
        # Integration tests
    end
end
```

### Test Documentation
Document test purpose and expected behavior:

```julia
@testset "Soliton propagation" begin
    # Test that a fundamental soliton maintains its shape
    # during propagation in a lossless fiber
    
    # Setup
    A₀ = sech_pulse(t, power, duration)
    wg = Waveguide(0.0, [-20e-27], 0.1, 1550e-9, 1.0)
    
    # Propagate
    sol = gnlse(A₀, t, wg)
    
    # Verify shape preservation
    @test isapprox(abs2.(sol.At[1,:]), abs2.(sol.At[end,:]), 
                   rtol=1e-3)
end
```
