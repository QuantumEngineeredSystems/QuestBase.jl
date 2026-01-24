using Test
using Symbolics
using SymbolicUtils: Fixpoint, Prewalk, PassThrough, BasicSymbolic

using QuestBase: @eqtest

@testset "exp(x)^n => exp(x*n)" begin
    using QuestBase: expand_all, expand_exp_power
    @variables a n

    @eqtest expand_exp_power(exp(a)^3) == exp(3 * a)
    @eqtest simplify(exp(a)^3) == exp(3 * a)
    @eqtest simplify(exp(a)^n) == exp(n * a)
    @eqtest expand_all(exp(a)^3) == exp(3 * a)
    @eqtest expand_all(exp(a)^3) == exp(3 * a)
    @eqtest expand_all(im * exp(a)^5) == im * exp(5 * a)
end

@testset "exp(a)*exp(b) => exp(a+b)" begin
    using QuestBase: simplify_exp_products
    @variables a b

    @eqtest simplify_exp_products(exp(a) * exp(b)) == exp(a + b)
    @eqtest simplify_exp_products(exp(3a) * exp(4b)) == exp(3a + 4b)
    @eqtest simplify_exp_products(im * exp(3a) * exp(4b)) == im * exp(3a + 4b)
end

@testset "euler" begin
    @variables a b
    @eqtest cos(a) + im * sin(a) == exp(im * a)
    @eqtest exp(a) * cos(b) + im * sin(b) * exp(a) == exp(a + im * b)
end

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
    @test string(eq.rhs) == "a" broken = true

    @eqtest drop_powers([a^2 + a + b, b], a, 2) == [a + b, b]
    @eqtest drop_powers([a^2 + a + b, b], [a, b], 2) == [a + b, b]
end

@testset "trig_to_exp and trig_to_exp" begin
    using QuestBase: expand_all, trig_to_exp, exp_to_trig
    @testset "Num" begin
        @variables f t
        cos_euler(x) = (exp(im * x) + exp(-im * x)) / 2
        sin_euler(x) = (exp(im * x) - exp(-im * x)) / (2 * im)

        # automatic conversion between trig and exp form
        trigs = [cos(f * t), sin(f * t)]
        for (i, trig) in pairs(trigs)
            z = trig_to_exp(trig)
            @eqtest simplify(expand(exp_to_trig(z))) == trig
        end
        trigs′ = [cos_euler(f * t), sin_euler(f * t)]
        for (i, trig) in pairs(trigs′)
            z = trig_to_exp(trig)
            @eqtest simplify(expand(exp_to_trig(z))) == trigs[i]
        end
    end

    # @testset "BasicSymbolic" begin
    #     @syms f t
    #     cos_euler(x) = (exp(im * x) + exp(-im * x)) / 2
    #     sin_euler(x) = (exp(im * x) - exp(-im * x)) / (2 * im)

    #     # automatic conversion between trig and exp form
    #     trigs = [cos(f * t), sin(f * t)]
    #     for (i, trig) in pairs(trigs)
    #         z = trig_to_exp(trig)
    #         @eqtest expand(exp_to_trig(z)) == trig
    #     end
    #     trigs′ = [cos_euler(f * t), sin_euler(f * t)]
    #     for (i, trig) in pairs(trigs′)
    #         z = trig_to_exp(trig)
    #         @eqtest expand(exp_to_trig(z)) == trigs[i]
    #     end
    # end
end

@testset "harmonic" begin
    using QuestBase: is_harmonic

    @variables a, b, c, t, x(t), f, y(t)

    @test is_harmonic(cos(f * t), t)
    @test is_harmonic(1, t)
    @test !is_harmonic(cos(f * t^2 + a), t)
    @test !is_harmonic(a + t, t)

    using QuestBase: DifferentialEquation

    dEOM = DifferentialEquation([a + x, t^2 + cos(t)], [x, y])
    @test !is_harmonic(dEOM, t)
end

@testset "fourier" begin
    using QuestBase: fourier_cos_term, fourier_sin_term
    using QuestBase.Symbolics: expand

    @variables f t θ a b

    @eqtest fourier_cos_term(cos(f * t)^2, f, t) == 0
    @eqtest fourier_sin_term(sin(f * t)^2, f, t) == 0

    @eqtest fourier_cos_term(cos(f * t)^2, 2 * f, t) == 1//2
    @eqtest fourier_sin_term(cos(f * t)^2, 2 * f, t) == 0
    @eqtest fourier_cos_term(sin(f * t)^2, 2 * f, t) == -1//2
    @eqtest fourier_sin_term(sin(f * t)^2, 2 * f, t) == 0

    @eqtest fourier_cos_term(cos(f * t), f, t) == 1
    @eqtest fourier_sin_term(sin(f * t), f, t) == 1

    @eqtest fourier_cos_term(cos(f * t + θ), f, t) == cos(θ)
    @eqtest fourier_sin_term(cos(f * t + θ), f, t) == -sin(θ)

    term =
        (a * sin(f * t) + b * cos(f * t)) *
        (a * sin(2 * f * t) + b * cos(2 * f * t)) *
        (a * sin(3 * f * t) + b * cos(3 * f * t))

    @eqtest fourier_cos_term(term, 2 * f, t) == expand(1//4 * (a^2 * b + b^3))
    @eqtest fourier_cos_term(term, 4 * f, t) == expand(1//4 * (a^2 * b + b^3))
    @eqtest fourier_cos_term(term, 6 * f, t) == expand(1//4 * (-3 * a^2 * b + b^3))
    @eqtest fourier_sin_term(term, 2 * f, t) == expand(1//4 * (a^3 + a * b^2))
    @eqtest fourier_sin_term(term, 4 * f, t) == expand(1//4 * (a^3 + a * b^2))
    @eqtest fourier_sin_term(term, 6 * f, t) == expand(1//4 * (-a^3 + 3 * a * b^2))

    # try something harder!
    term = (a + b * cos(f * t + θ)^2)^3 * sin(f * t)
    @eqtest fourier_sin_term(term, f, t) == expand(
        a^3 + a^2 * b * 3//2 + 9//8 * a * b^2 + 5//16 * b^3 -
        3//64 * b * (16 * a^2 + 16 * a * b + 5 * b^2) * cos(2 * θ),
    )

    @eqtest fourier_cos_term(term, f, t) ==
        expand(-3//64 * b * (16 * a^2 + 16 * a * b + 5 * b^2) * sin(2 * θ))

    # FTing at zero : picking out constant terms
    @eqtest fourier_cos_term(cos(f * t), 0, t) == 0
    @eqtest fourier_cos_term(cos(f * t)^3 + 1, 0, t) == 1
    @eqtest fourier_cos_term(cos(f * t)^2 + 1, 0, t) == 3//2
    @eqtest fourier_cos_term((cos(f * t)^2 + cos(f * t))^3, 0, t) == 23//16
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
        @test simplify_complex(z).val isa BasicSymbolic
    end

    z = Complex{Num}(1 + 0 * im)
    @test SymbolicUtils.is_literal_number(simplify_complex(z).val)
end

@testset "get_all_terms" begin
    using QuestBase: get_all_terms
    @variables a, b, c

    @eqtest sort(get_all_terms(a + b + c); by=string) == sort([a, b, c]; by=string)
    @eqtest sort(get_all_terms(a * b * c); by=string) == sort([a, b, c]; by=string)
    @eqtest sort(get_all_terms(a / b); by=string) == sort([a, b]; by=string)
    @eqtest sort(get_all_terms(a^2 + b^2 + c^2); by=string) == sort([a^2, b^2, c^2]; by=string)
    @eqtest sort(get_all_terms(a^2 / b^2); by=string) == sort([a^2, b^2]; by=string)
    @eqtest sort(get_all_terms(2 * b^2); by=string) == sort([2, b^2]; by=string)
    @eqtest sort(get_all_terms(2 * b^2 ~ a); by=string) == sort([2, b^2, a]; by=string)
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
    @test string.(expand_fraction(d * (a + b) / c)) == "d*a /c + d*b / c + d / c" broken =
        true
end
@testset "count_derivatives" begin
    using QuestBase: count_derivatives, d
    @variables t x(t) y(t)
    @test count_derivatives(x) == 0
    @test count_derivatives(d(x, t)) == 1
    @test count_derivatives(d(d(x, t), t)) == 2
    @test_throws ErrorException count_derivatives(d(d(5 * x, t), t))
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
