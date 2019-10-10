using Documenter, Tweedie

makedocs(;
    modules=[Tweedie],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jkbest2/Tweedie.jl/blob/{commit}{path}#L{line}",
    sitename="Tweedie.jl",
    authors="John Best",
    assets=String[],
)

deploydocs(;
    repo="github.com/jkbest2/Tweedie.jl",
)
