CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using QuestBase
using Documenter

include("pages.jl")

makedocs(;
    sitename="QuestBase.jl",
    authors="Quest group",
    modules=QuestBase,
    format=Documenter.HTML(;
        canonical="https://quantumengineeredsystems.github.io/QuestBase.jl/stable/"
    ),
    pages=pages,
    clean=true,
    linkcheck=true,
    warnonly=:missing_docs,
    draft=!CI,
    doctest=false,  # We test it in the CI, no need to run it here
)

if CI
    deploydocs(;
        repo="github.com/QuantumEngineeredSystems/QuestBase.jl";
        # devbranch="master",
        # target="build",
        # branch="gh-pages",
        push_preview=true,
    )
end
