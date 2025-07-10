# Contributing

We welcome contributions to FiberNlse! Whether you want to report a bug, suggest a feature, or contribute code, your input is valuable.

## How to Contribute

1. **Report Issues**: If you encounter any issues, please open an issue on our [GitHub repository](https://github.com/brian-sinquin/FiberNlse.jl/issues).
2. **Suggest Features**: Have an idea for a new feature? Open a feature request issue.
3. **Contribute Code**:
   - Fork the repository.
   - Create a new branch for your feature or bugfix.
   - Submit a pull request with a clear description of your changes.

## Development Setup

To set up a development environment:

```julia
using Pkg
Pkg.develop(url="https://github.com/brian-sinquin/FiberNlse.jl")
Pkg.instantiate()
```

## Code Style

We follow the [BlueStyle](https://github.com/invenia/BlueStyle) coding conventions. Please ensure your code adheres to these guidelines.

## Testing

Before submitting a pull request, ensure that all tests pass:

```julia
using Pkg
Pkg.test("FiberNlse")
```

Thank you for contributing to FiberNlse!
