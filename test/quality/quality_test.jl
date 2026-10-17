using QuestBase
using Test
using Aqua
using CheckConcreteStructs: all_concrete
using ExplicitImports

@testset "Concretely typed" begin
    all_concrete(QuestBase.HarmonicVariable)
    all_concrete(QuestBase.HarmonicEquation)
    all_concrete(QuestBase.DifferentialEquation)
end

@testset "Code quality" begin
    @test check_no_stale_explicit_imports(QuestBase) == nothing
    @test check_all_explicit_imports_via_owners(QuestBase) == nothing
    Aqua.test_ambiguities([QuestBase])
    Aqua.test_all(QuestBase;)
end
