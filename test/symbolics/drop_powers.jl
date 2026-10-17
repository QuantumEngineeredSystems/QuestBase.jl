using Test
using Symbolics
using SymbolicUtils: BasicSymbolic
using QuestBase: QuestBase, @eqtest

@testset "powers" begin
    using QuestBase: drop_powers, max_power
    using QuestBase.Symbolics: expand

    @variables a, b, c

    @eqtest max_power(a^2 + b, a) == 2
    @eqtest max_power(a * ((a + b)^4)^2 + a, a) == 9
    @eqtest max_power([a * ((a + b)^4)^2 + a, a^2], a) == 9
    @eqtest max_power(a + im * a^2, a) == 2

    @eqtest drop_powers(a^2 + b, a, 1) == b
    @eqtest drop_powers((a + b)^2, a, 1) == b^2
    @eqtest drop_powers((a + b)^2, [a, b], 1) == 0
    @eqtest drop_powers((a + b)^3 + (a + b)^5, [a, b], 4) == expand((a + b)^3)
    a^2 + a ~ a
    eq = drop_powers(a^2 + a ~ a, a, 2)
    @eqtest [eq.lhs, eq.rhs] == [a, a]
    # eq = drop_powers(a^2 + a ~ b, [a, b], 2) # broken
    @eqtest [eq.lhs, eq.rhs] == [a, a]
    eq = drop_powers(a^2 + a + b ~ a, a, 2)
    @test string(eq.rhs) == "a"

    @eqtest drop_powers([a^2 + a + b, b], a, 2) == [a + b, b]
    @eqtest drop_powers([a^2 + a + b, b], [a, b], 2) == [a + b, b]
end

@testset "powers of raw symbolic variables" begin
    using QuestBase: drop_powers
    using SymbolicUtils: unwrap
    @variables a b
    raw_a, raw_b = unwrap(a), unwrap(b)
    expr = unwrap((a + b)^3)
    expected = expand(3a^2 * b + 3a * b^2 + b^3)

    result = drop_powers(expr, raw_a, 3)
    @test result isa BasicSymbolic
    @eqtest result == expected
    @eqtest drop_powers(expr, [raw_a], 3) == expected
    @test isequal(drop_powers(expr, [raw_a, raw_b], 3), unwrap(Num(0)))
    @eqtest drop_powers([expr, raw_b], raw_a, 3) == [expected, b]
    @eqtest drop_powers(a^3 + b, raw_a, 3) == b
    @eqtest drop_powers(expr, a, 3) == expected
    @test isequal(drop_powers(expr, [a, raw_b], 3), unwrap(Num(0)))

    eq = drop_powers(a^3 + a ~ b, raw_a, 3)
    @eqtest [eq.lhs, eq.rhs] == [a, b]
    @test_throws ArgumentError drop_powers(expr, [1], 3)
    @test_throws ArgumentError drop_powers(expr, raw_a, -1)
    @test isequal(drop_powers(expr, raw_a, 0), unwrap(Num(0)))
    @eqtest drop_powers(Num(7), a, 0) == 0

    @variables epsilon ϵ
    @eqtest drop_powers(epsilon + ϵ^2 + a^3, a, 3) == epsilon + ϵ^2
end
