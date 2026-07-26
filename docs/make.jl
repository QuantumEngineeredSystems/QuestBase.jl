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
    # HarmonicBalance.jl's docs site currently 404s below its root, so a broken external
    # link must not fail the build
    warnonly=[:missing_docs, :linkcheck],
    draft=(!CI),
    doctest=false,  # We test it in the CI, no need to run it here
)

if CI
    deploydocs(;
        repo="github.com/QuantumEngineeredSystems/QuestBase.jl",
        devbranch="main",
        target="build",
        branch="gh-pages",
        push_preview=true,
    )
end
