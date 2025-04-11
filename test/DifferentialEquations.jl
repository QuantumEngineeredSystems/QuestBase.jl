
using Test
using Symbolics
using OrderedCollections
using QuestBase:
    d,
    DifferentialEquation,
    @eqtest,
    get_variables,
    is_harmonic,
    get_independent_variables,
    add_harmonic!,
    rearrange_standard!,
    is_rearranged_standard,
    rearrange,
    substitute_all,
    _remove_brackets,
    declare_variables,
    dummy_symbolic_Jacobian,
    get_all_terms

# Setup common variables
@variables t x(t) y(t) ω0 ω F k

@testset "Constructor Tests" begin
    # Single equation constructor
    diff_eq1 = DifferentialEquation(d(x, t, 2) + ω0^2 * x ~ F * cos(ω * t), x)
    @test length(diff_eq1.equations) == 1
    @eqtest first(keys(diff_eq1.equations)) == x

    # Multiple equations constructor
    eqs =
        [d(x, t, 2) + ω0^2 * x - k * y, d(y, t, 2) + ω0^2 * y - k * x] .~
        [F * cos(ω * t), 0]
    vars = [x, y]
    diff_eq2 = DifferentialEquation(eqs, vars)
    @test length(diff_eq2.equations) == 2
    @eqtest collect(keys(diff_eq2.equations)) == vars

    # Expression to equation constructor
    expr = d(x, t, 2) + ω0^2 * x
    diff_eq3 = DifferentialEquation(expr, x)
    @test length(diff_eq3.equations) == 1
    @test diff_eq3.equations[x].rhs == 0
end

@testset "Helper Functions" begin
    diff_eq = DifferentialEquation(d(x, t, 2) + ω0^2 * x ~ F * cos(ω * t), x)

    # Test get_variables
    @eqtest get_variables(diff_eq) == [x]

    # Test is_harmonic
    @test is_harmonic(diff_eq, t)

    # Test get_independent_variables
    @eqtest [t] == get_independent_variables(diff_eq)
end

@testset "Harmonic Manipulation" begin
    diff_eq = DifferentialEquation(d(x, t, 2) + ω0^2 * x ~ F * cos(ω * t), x)

    # Test add_harmonic!
    @test isempty(diff_eq.harmonics[x])
    add_harmonic!(diff_eq, x, ω)
    @eqtest ω == first(diff_eq.harmonics[x])
end

@testset "Equation Rearrangement" begin
    diff_eq = DifferentialEquation(d(x, t, 2) + ω0^2 * x ~ F * cos(ω * t), x)

    # Test is_rearranged_standard
    @test !is_rearranged_standard(diff_eq)

    # Test rearrange_standard!
    rearrange_standard!(diff_eq)
    @test is_rearranged_standard(diff_eq)

    # Test rearrange
    new_eq = rearrange(diff_eq, [x])
    @test new_eq !== diff_eq  # Should be a new instance
    @eqtest new_eq != diff_eq
end

@testset "Error Cases" begin
    @testset "Equation{Vector}" begin
        # define equation of motion
        @variables ω1, ω2, t, ω, F, γ, α1, α2, k, x(t), y(t)
        rhs = [
            d(x, t, 2) + ω1^2 * x + γ * d(x, t) + α1 * x^3 - k * y,
            d(d(y, t), t) + ω2^2 * y + γ * d(y, t) + α2 * y^3 - k * x,
        ]
        eqs = rhs .~ [F * cos(ω * t), 0]

        @test eqs != (rhs ~ [F * cos(ω * t), 0])
        @test eqs == (rhs .~ [F * cos(ω * t), 0])
        @test_throws ArgumentError DifferentialEquation(rhs ~ [F * cos(ω * t), 0], [x, y])
        @test_throws ArgumentError DifferentialEquation([rhs[1]] ~ [F * cos(ω * t)], [x, y])
    end

    # Test invalid equation format
    @test_throws ArgumentError DifferentialEquation([d(x, t, 2) ~ F * cos(ω * t)], [x, y])
end
