@testset "Concretely typed" begin
    using QuestBase

    using CheckConcreteStructs

    all_concrete(QuestBase.HarmonicVariable)
    all_concrete(QuestBase.HarmonicEquation)
    all_concrete(QuestBase.DifferentialEquation)
    all_concrete(QuestBase.HarmonicVariable)
end

@testset "Code linting" begin
    using JET
    rep = report_package(QuestBase; target_modules=(QuestBase,))
    @show rep
    @test length(JET.get_reports(rep)) == 0
end

@testset "Code quality" begin
    using QuestBase
    using ExplicitImports, Aqua

    @test check_no_stale_explicit_imports(QuestBase) == nothing
    @test check_all_explicit_imports_via_owners(QuestBase) == nothing
    Aqua.test_ambiguities([QuestBase])
    Aqua.test_all(QuestBase;)
end
