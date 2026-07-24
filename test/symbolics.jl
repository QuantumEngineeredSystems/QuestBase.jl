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
    # In SymbolicUtils v4, exp(im*x) stays as a symbolic Term and does NOT
    # auto-decompose to cos(x)+im*sin(x). Test via QuestBase's exp_to_trig.
    using QuestBase: exp_to_trig
    using SymbolicUtils: @syms
    @syms a_eu
    @eqtest expand(exp_to_trig(exp(im * a_eu))) == cos(a_eu) + im * sin(a_eu)
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
    @test string(eq.rhs) == "a"

    @eqtest drop_powers([a^2 + a + b, b], a, 2) == [a + b, b]
    @eqtest drop_powers([a^2 + a + b, b], [a, b], 2) == [a + b, b]
end

@testset "trig_to_exp and exp_to_trig" begin
    using QuestBase: expand_all, trig_to_exp, exp_to_trig
    @testset "Num" begin
        @variables f t

        # automatic conversion between trig and exp form
        trigs = [cos(f * t), sin(f * t)]
        for (i, trig) in pairs(trigs)
            z = trig_to_exp(trig)
            @eqtest expand(exp_to_trig(z)) == trig
        end
    end

    @testset "BasicSymbolic" begin
        using SymbolicUtils: @syms
        @syms f_bs t_bs

        # At BasicSymbolic level, exp(im*x) stays symbolic - roundtrip works cleanly
        trigs = [cos(f_bs * t_bs), sin(f_bs * t_bs)]
        for (i, trig) in pairs(trigs)
            z = trig_to_exp(trig)
            @eqtest expand(exp_to_trig(z)) == trig
        end
    end
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

@testset "trig_reduce" begin
    using QuestBase: trig_reduce

    @variables f t a ω

    # the Pythagorean identity collapses in a bare expression …
    @eqtest trig_reduce(cos(f * t)^2 + sin(f * t)^2) == 1

    # … and, crucially, in the denominator of a combined fraction. The exponential
    # round-trip only touches the numerator, so without `reduce_denominator` the
    # identity survives in the denominator, `get_independent` sees a time-dependent
    # denominator and discards the whole fraction, collapsing the Krylov-Bogoliubov
    # slow-flow equations to `0 ~ d/dT`.
    @eqtest trig_reduce(a / (cos(f * t)^2 + sin(f * t)^2)) == a
    @eqtest trig_reduce(a / (ω * cos(f * t)^2 + ω * sin(f * t)^2)) == a / ω

    # trig-free denominators pass through `reduce_denominator` untouched
    using QuestBase: reduce_denominator
    @eqtest reduce_denominator(a / (ω^2 + 1)) == a / (ω^2 + 1)
end

@testset "rearrange! reduces determinant denominators" begin
    # `symbolic_linear_solve` returns nested fractions whose denominators hold the
    # coefficient determinant; for a trigonometric ansatz it contains cos² + sin² = 1.
    # `rearrange!` must flatten and collapse it, otherwise every consumer works with
    # superficially time-dependent denominators (this stalled the Krylov-Bogoliubov
    # order-2 averaging for hours).
    using QuestBase: HarmonicEquation, HarmonicVariable, DifferentialEquation
    using QuestBase: rearrange!, d
    @variables ω t T u1(T) v1(T) x(t)

    eqs = [
        cos(ω * t) * d(u1, T) + sin(ω * t) * d(v1, T) ~ u1,
        -sin(ω * t) * d(u1, T) + cos(ω * t) * d(v1, T) ~ v1,
    ]
    hvars = [HarmonicVariable(u1, "u", "u", ω, x), HarmonicVariable(v1, "v", "v", ω, x)]
    natural_eq = DifferentialEquation(d(x, t, 2) + x ~ 0, x)
    eom = HarmonicEquation(eqs, hvars, natural_eq)

    rearrange!(eom, d([u1, v1], T))
    lhs = QuestBase.Num.(getfield.(eom.equations, :lhs))
    @eqtest lhs[1] == u1 * cos(ω * t) - v1 * sin(ω * t)
    @eqtest lhs[2] == u1 * sin(ω * t) + v1 * cos(ω * t)
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
