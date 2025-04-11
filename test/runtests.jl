using QuestBase
using Test

using Random
const SEED = 0x8f88209c
Random.seed!(SEED)

@testset "Code quality" begin
    include("code_quality.jl")
end

@testset "Symbolics customised" begin
    include("symbolics.jl")
end

@testset "DifferentialEquations" begin
    include("DifferentialEquations.jl")
end

@testset "HarmonicVariable" begin
    include("HarmonicVariable.jl")
end

@testset "HarmonicEquation" begin
    include("HarmonicEquation.jl")
end
