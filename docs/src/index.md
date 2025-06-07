# FiberNlse

Documentation for [FiberNlse](https://github.com/brian-sinquin/FiberNlse.jl).

## Overview

FiberNlse is a Julia package for simulating pulse propagation in optical fibers using the generalized nonlinear Schrödinger equation (GNLSE).

## Documentation

- [Installation](installation.md)
- [Theoretical background](theory.md)
- [Usage](usage.md)
- [API Reference](api.md)
- [Examples](examples.md)

## Index

```@index
```

using Documenter, Literate

# Generate example notebooks/scripts
Literate.markdown("examples/propagation.jl", "src/examples")
Literate.markdown("examples/supercontinuum.jl", "src/examples")

makedocs(
    ...
    modules = [FiberNlse],
    ...
)