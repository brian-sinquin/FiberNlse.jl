# Installation

FiberNlse can be installed using Julia's package manager.

## Stable Release

Open Julia and run:

```julia
using Pkg
Pkg.add("FiberNlse")
```

## Development Version

To use the latest development version directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/brian-sinquin/FiberNlse.jl")
```

## Development Mode

If you want to contribute or make local changes, you can develop the package:

```julia
using Pkg
Pkg.develop(url="https://github.com/brian-sinquin/FiberNlse.jl")
```

Or, clone the repository and add it as a dev package:

```sh
git clone https://github.com/brian-sinquin/FiberNlse.jl.git
```
```julia
using Pkg
Pkg.develop(path="path/to/FiberNlse.jl")
```
