using Test
using Symbolics
using SymbolicUtils: BasicSymbolic
using QuestBase: QuestBase, @eqtest

@testset "harmonic" begin
    using QuestBase: is_harmonic, d

    @variables a, b, c, t, x(t), f, y(t)

    @test is_harmonic(cos(f * t), t)
    @test is_harmonic(1, t)
    @test !is_harmonic(cos(f * t^2 + a), t)
    @test !is_harmonic(a + t, t)

    # only the summands holding `t` are expanded; the rest must not change the verdict
    @test is_harmonic(a * x^3 + b * d(x, t, 2) + f * cos(f * t), t)
    @test is_harmonic(a * x + (b + cos(f * t))^2, t)
    @test is_harmonic(a * cos(f * t)^2 * x + sin(2f * t)^3, t)
    @test is_harmonic(cos(f * t)^(-1), t) # 1/cos: a denominator is split, the projection rejects it
    @test !is_harmonic(cos(f * t)^(1//2), t)
    @test !is_harmonic(cos(f * t^2)^2, t)
    @test !is_harmonic(a * x^3 + t * cos(f * t), t)
    @test is_harmonic(x + a / (b + cos(f * t)), t) # a denominator is split into its terms
    @test !is_harmonic(x^2 + cos(f * t)^2 * exp(t), t)

    using QuestBase: DifferentialEquation

    dEOM = DifferentialEquation([a + x, t^2 + cos(t)], [x, y])
    @test !is_harmonic(dEOM, t)
end

@testset "polynomial in the variables" begin
    using QuestBase: is_polynomial, nonpolynomial_terms, d, DifferentialEquation

    @variables a, b, t, x(t), y(t), ω

    @test is_polynomial(a * x^3 + b * d(x, t) * x^2 - a, [x])
    @test is_polynomial(a * cos(ω * t) * x - b * x * y, [x, y])
    @test is_polynomial(d(x, t, 2) + a * x / b, [x]) # variable-free denominator
    @test is_polynomial(a + cos(ω * t), [x]) # no variable at all
    @test is_polynomial(a * x^2.0, [x]) # integer-valued exponent of another type
    @test all(n -> is_polynomial(a * x^n + b * d(x, t) * x^n, [x]), 0:12)

    @test !is_polynomial(a * sin(x), [x])
    @test !is_polynomial(a * exp(x + y), [x, y])
    @test !is_polynomial(a / x, [x])
    @test !is_polynomial(x^(-2), [x])
    @test !is_polynomial(sqrt(x), [x])
    @test !is_polynomial(a^x, [x])
    @test !is_polynomial(d(sin(x), t), [x]) # hidden inside a derivative

    # the offending terms are reported, the harmless ones are not
    terms = nonpolynomial_terms(a * x^3 + a * sin(x) ~ b * cos(ω * t), [x])
    @eqtest Num.(terms) == [sin(x)]

    # a variable of another equation does not make a term non-polynomial
    @test is_polynomial(sin(y) + x, [x])

    dEOM = DifferentialEquation(
        [d(x, t, 2) + x + a * x^3, d(y, t, 2) + y + a * sin(y)], [x, y]
    )
    @eqtest Num.(nonpolynomial_terms(dEOM)) == [sin(y)]
    @test isempty(nonpolynomial_terms(DifferentialEquation([x + a * x^3, y], [x, y])))
end

@testset "_apply_termwise" begin
    using QuestBase: _apply_termwise

    @variables a, b, c

    @eqtest _apply_termwise(x -> x^2, a + b + c) == a^2 + b^2 + c^2
    @eqtest _apply_termwise(x -> x^2, a * b * c) == a^2 * b^2 * c^2
    @eqtest _apply_termwise(x -> x^2, a / b) == a^2 / b^2
end

@testset "simplify_complex" begin
    using QuestBase: simplify_complex
    @variables a, b, c
    for z in Complex{Num}[a, a * b, a / b]
        @test Symbolics.unwrap(simplify_complex(z)) isa BasicSymbolic
    end

    z = Complex{Num}(1 + 0 * im)
    @test simplify_complex(z) isa Number
end

@testset "get_all_terms" begin
    using QuestBase: get_all_terms
    @variables a, b, c

    # Use Set comparison since ordering may vary across SymbolicUtils versions
    @test Set(get_all_terms(a + b + c)) == Set([a, b, c])
    @test Set(get_all_terms(a * b * c)) == Set([a, b, c])
    @test Set(get_all_terms(a / b)) == Set([a, b])
    @test Set(get_all_terms(a^2 + b^2 + c^2)) == Set([a^2, b^2, c^2])
    @test Set(get_all_terms(a^2 / b^2)) == Set([a^2, b^2])
    @test Set(get_all_terms(2 * b^2)) == Set([2, b^2])
    @test Set(get_all_terms(2 * b^2 ~ a)) == Set([2, b^2, a])
end

@testset "get_independent" begin
    using QuestBase: get_independent
    @variables a, b, c, t

    @test get_independent(a + b + c, t) isa Num
    @eqtest get_independent(a + b + c, t) == a + b + c
    @eqtest get_independent(a * b * c, t) == a * b * c
    @eqtest get_independent(a / b, t) == a / b
    @eqtest get_independent(a^2 + b^2 + c^2, t) == a^2 + b^2 + c^2
    @eqtest get_independent(a^2 / b^2, t) == a^2 / b^2
    @eqtest get_independent(2 * b^2, t) == 2 * b^2
    @eqtest get_independent(cos(t), t) == 0
    @eqtest get_independent(cos(t)^2 + 5, t) == 5
    @eqtest get_independent(a + im * b, t) == a + im * b
    @eqtest get_independent(a + im * b + cos(t)^2 + 5, t) == 5 + a + im * b
end

@testset "expand_fraction" begin
    using QuestBase: expand_fraction
    @variables a, b, c, d

    @eqtest expand_fraction((a + b) / c) == a / c + b / c
    @eqtest expand_fraction(d * (a + b) / c) == a * d / c + b * d / c
end
@testset "count_derivatives" begin
    using QuestBase: count_derivatives, d
    @variables t x(t) y(t)
    @test count_derivatives(x) == 0
    @test count_derivatives(d(x, t)) == 1
    @test count_derivatives(d(d(x, t), t)) == 2
    # In Symbolics v7, Differential stores order directly, so d(d(5*x,t),t)
    # becomes Differential(t,2)(5*x(t)) which is still a valid derivative term
    @test count_derivatives(d(d(5 * x, t), t)) == 2
end

@testset "substitute_all" begin
    using QuestBase: substitute_all
    @variables a, b, c, d, e, f, g, h

    pairs = (a => b, c => d, e => f, g => h)
    rules = Dict(a => b, c => d, e => f, g => h)

    @eqtest substitute_all(a + b + c + d + e + f + g + h, rules) ==
        2 * b + 2 * d + 2 * f + 2 * h
    @eqtest substitute_all(a + b, a => b) == 2 * b
    @eqtest substitute_all(a * b * c * d * e * f * g * h, rules) == b^2 * d^2 * f^2 * h^2
    @eqtest substitute_all([a, c, e], rules) == [b, d, f]
    @eqtest substitute_all(a + b * im, rules) == b + b * im
end
