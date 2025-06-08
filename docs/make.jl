using Documenter, Literate

# Generate example notebooks/scripts
Literate.markdown("examples/supercontinuum_example.jl", "src/examples"; documenter=true)
Literate.markdown("examples/soliton_example.jl", "src/examples"; documenter=true)

makedocs(;
    sitename="FiberNlse.jl",
    modules=[FiberNlse],
    checkdocs=:exports,
    authors="curio-sitas <brian.sinquin@gmail.com> and contributors",
    format=Documenter.HTML(;
        canonical="https://brian-sinquin.github.io/FiberNlse.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "Supercontinuum" => "examples/supercontinuum_example.md",
            "Soliton" => "examples/soliton_example.md",
        ],
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/brian-sinquin/FiberNlse.jl",
    push_preview=true,
)
