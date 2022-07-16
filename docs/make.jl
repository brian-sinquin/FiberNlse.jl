using Documenter
using FiberNlse

makedocs(
    sitename = "FiberNlse",
    format = Documenter.HTML(),
    modules = [FiberNlse]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/curio-sitas/FiberNlse.jl.git"
)
