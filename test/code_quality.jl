@testset "Concretely typed" begin
    using QuestBase

    using CheckConcreteStructs

    all_concrete(QuestBase.HarmonicVariable)
    all_concrete(QuestBase.HarmonicEquation)
    all_concrete(QuestBase.DifferentialEquation)
    all_concrete(QuestBase.HarmonicVariable)
end

CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
@testset "Code linting" begin
    # JET is skipped on CI because Moshi.@match generates complex pattern-matching
    # dispatch code that causes JET analysis to time out.
    if !CI
        using JET
        rep = report_package(QuestBase; target_modules=(QuestBase,))
        @show rep
        # Filter out Moshi @match false positives: generated variable bindings
        # in pattern match arms trigger "may be undefined" warnings in JET
        real_reports = filter(JET.get_reports(rep)) do r
            !contains(string(r), "##Call#") && !contains(string(r), "##And#")
        end
        @test length(real_reports) == 0
    end
end

@testset "Code quality" begin
    using QuestBase
    using ExplicitImports, Aqua

    @test check_no_stale_explicit_imports(QuestBase) == nothing
    @test check_all_explicit_imports_via_owners(QuestBase) == nothing
    Aqua.test_ambiguities([QuestBase])
    Aqua.test_all(QuestBase;)
end
