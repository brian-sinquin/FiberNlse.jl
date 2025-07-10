# Add the src directory to LOAD_PATH
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using Documenter, FiberNlse

cd(@__DIR__)

makedocs(;
    sitename="FiberNlse.jl",
    modules=[FiberNlse],
    checkdocs=:exports,
    authors="brian-sinquin <brian.sinquin@gmail.com> and contributors",
    format=Documenter.HTML(;
        canonical="https://brian-sinquin.github.io/FiberNlse.jl",
        edit_link="main",
        prettyurls=get(ENV, "CI", nothing) == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => [
            "Installation" => "installation.md",
            "Quick Start Guide" => "quickstart.md",
        ],
        "User Guide" => [
            "Usage" => "usage.md",
            "Theoretical Background" => "theory.md",
            "Visualization Guide" => "visualization.md",
        ],
        "Examples" => [
            "Basic Effects" => [
                "Soliton Propagation" => "examples/soliton.md",
                "Soliton Fission" => "examples/soliton_fission.md",
                "Pulse Compression" => "examples/compression.md",
            ],
            "Advanced Effects" => [
                "Supercontinuum Generation" => "examples/supercontinuum.md",
            ],
            "Different Waveguides" => "examples/waveguides.md",
        ],
        "Technical Reference" => "reference.md",
        "API Reference" => "api.md",
        "Developer Guide" => [
            "Contributing" => "developer/contributing.md",
            "Code Style" => "developer/style.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/brian-sinquin/FiberNlse.jl",
    push_preview=true,
)
