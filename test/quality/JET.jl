using QuestBase
using Test
using JET

@testset "JET" begin
    JET.test_package(QuestBase; target_modules=(QuestBase,))
end
