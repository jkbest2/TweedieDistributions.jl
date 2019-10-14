using Documenter, Tweedie

makedocs(;
    modules=[TweedieDistributions],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jkbest2/TweedieDistributions.jl/blob/{commit}{path}#L{line}",
    sitename="TweedieDistributions.jl",
    authors="John Best",
    assets=String[],
)

deploydocs(;
    repo="github.com/jkbest2/TweedieDistributions.jl",
)
