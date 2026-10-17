using Test
using Symbolics
using SymbolicUtils: Fixpoint, Prewalk, PassThrough, BasicSymbolic

using QuestBase: QuestBase, @eqtest

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

    # a product stores a repeated factor as a power, so `exp(a)*exp(a)` never reaches the
    # merge as two separate `exp` factors. It arrives as the single factor `exp(a)^2`.
    @eqtest simplify_exp_products(exp(a) * exp(a)) == exp(2a)
    @eqtest simplify_exp_products(b * exp(a) * exp(a)) == b * exp(2a)
    @eqtest simplify_exp_products(exp(a)^2 * exp(b)) == exp(2a + b)
    # and the exponents must cancel to nothing, as they do for two plain factors
    @eqtest simplify_exp_products(exp(a)^2 * exp(-2a) * b) == b
    @eqtest simplify_exp_products(exp(a)^2 * exp(-a)^2) == 1
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

@testset "batched Fourier projection with phase shifts" begin
    using QuestBase: fourier_terms, fourier_cos_term, fourier_sin_term
    @variables ω t θ ψ F λ u v
    specs = [(ω, cos), (ω, sin), (3 * ω, cos), (3 * ω, sin), (0, cos)]
    drive = F * cos(ω * t + θ)
    pump = λ * (u * cos(ω * t) + v * sin(ω * t)) * cos(2 * ω * t + ψ)

    drive_coefficients = fourier_terms(drive, t, specs)
    @eqtest drive_coefficients == [F * cos(θ), -F * sin(θ), 0, 0, 0]
    pump_coefficients = fourier_terms(pump, t, specs)
    @eqtest pump_coefficients[1] == expand(λ * (u * cos(ψ) - v * sin(ψ)) / 2)
    @eqtest pump_coefficients[2] == expand(-λ * (u * sin(ψ) + v * cos(ψ)) / 2)

    expression = drive + pump + 2
    reference = [
        if f === cos
            fourier_cos_term(expression, frequency, t)
        else
            fourier_sin_term(expression, frequency, t)
        end for (frequency, f) in specs
    ]
    @eqtest fourier_terms(expression, t, specs) == reference
    equations = fourier_terms(expression ~ drive, t, specs)
    @eqtest getfield.(equations, :lhs) == reference
    @eqtest getfield.(equations, :rhs) == drive_coefficients
end

@testset "Fourier projection matches the discrete Fourier sum" begin
    using QuestBase: fourier_terms
    @variables ω t θ ψ F λ u v a b c m
    specs = [(ω, cos), (ω, sin), (2ω, cos), (2ω, sin), (3ω, cos), (3ω, sin), (0, cos)]
    harmonic(k) = a * cos(k * ω * t) + b * sin(k * ω * t)
    expressions = [
        harmonic(1)^3,
        harmonic(1)^5 + c * harmonic(1)^2,
        harmonic(1) * harmonic(2) * harmonic(3),
        0.3 * harmonic(1)^3 + 1.7 * harmonic(2),
        sin(ω * t)^4 * cos(3ω * t) + 7,
        F * cos(ω * t + θ),
        λ * (u * cos(ω * t) + v * sin(ω * t)) * cos(2ω * t + ψ),
        (a + b * cos(ω * t + θ)^2)^3 * sin(ω * t),
        cos(ω * t + θ)^3 * sin(2ω * t + ψ),
        cos(ω * t + θ - ψ)^2 * sin(ω * t + ψ),
        cos(ω * t + 1//3) * harmonic(1)^2,
        F * cos(ω * t + θ + 1//3) * harmonic(1),
        sin(ω * t - θ + 2)^2 * cos(ω * t),
        cos(ψ)^2 * cos(ω * t),
        (c * harmonic(1)^3 + harmonic(2)) / m,
        a / (b + c) * cos(ω * t)^3 + harmonic(1) / m,
        Num(5),
        Num(0),
    ]
    values = Dict(
        ω => 1.0,
        θ => 0.7,
        ψ => -0.4,
        F => 1.3,
        λ => 0.9,
        u => 0.2,
        v => -1.1,
        a => 0.8,
        b => -0.6,
        c => 1.7,
        m => 2.3,
    )
    numeric(x, point) = Symbolics.value(substitute(x, point; fold=Val(true)))
    # The rectangle rule over one period is exact for trig polynomials of degree below N.
    N = 64
    times = 2π .* (0:(N - 1)) ./ N
    for expression in expressions
        samples = [numeric(expression, merge(values, Dict(t => τ))) for τ in times]
        coefficients = fourier_terms(expression, t, specs)
        for ((k, f), coefficient) in zip(specs, coefficients)
            weight = isequal(k, 0) ? 1 / N : 2 / N
            expected = weight * sum(samples .* f.(numeric(k, values) .* times))
            @test numeric(coefficient, values) ≈ expected atol = 1e-12
        end
    end
end

@testset "Fourier projection is exact" begin
    using QuestBase: fourier_terms
    @variables ω ν t θ a b m
    specs = [(ω, cos), (ω, sin), (3ω, cos), (3ω, sin), (0, cos)]
    cubic = (a * cos(ω * t) + b * sin(ω * t))^3
    @eqtest fourier_terms(cubic, t, specs) == [
        3//4 * a^3 + 3//4 * a * b^2,
        3//4 * a^2 * b + 3//4 * b^3,
        1//4 * a^3 - 3//4 * a * b^2,
        3//4 * a^2 * b - 1//4 * b^3,
        0,
    ]

    # Phases enter through cos and sin of the phase, with the sign of the phase fixed and
    # numeric offsets kept exact.
    @eqtest fourier_terms(cos(ω * t + θ)^2, t, [(2ω, cos), (2ω, sin), (0, cos)]) ==
        [1//2 * cos(2θ), -1//2 * sin(2θ), 1//2]
    @eqtest fourier_terms(sin(ω * t - θ + 1//3), t, [(ω, cos)]) == [-sin(θ - 1//3)]

    # A t-free denominator divides each term.
    @eqtest fourier_terms(a * cos(ω * t)^3 / m, t, [(ω, cos)]) == [3//4 * a / m]

    # Incommensurate tones: products land on ω ± ν, never on ω alone.
    two_tone = (a * cos(ω * t) + b * cos(ν * t))^3
    tone_specs = [(ω, cos), (ν, cos), (ω + ν, cos), (2ω - ν, cos)]
    @eqtest fourier_terms(two_tone, t, tone_specs) ==
        [3//4 * a^3 + 3//2 * a * b^2, 3//2 * a^2 * b + 3//4 * b^3, 0, 3//4 * a^2 * b]

    equations = fourier_terms(a * cos(ω * t) ~ b * cos(ω * t)^3, t, [(ω, cos)])
    @eqtest only(equations).lhs == a
    @eqtest only(equations).rhs == 3//4 * b
end

@testset "Fourier projection rejects inputs without a finite Fourier series" begin
    using QuestBase: fourier_terms
    @variables ω t a
    @variables x(t)
    specs = [(ω, cos)]
    for expression in [
        a / (2 + cos(ω * t)),
        t * cos(ω * t),
        cos(ω * t^2),
        sqrt(cos(ω * t)),
        exp(cos(ω * t)),
        x * cos(ω * t),
        im * cos(ω * t),
    ]
        @test_throws ArgumentError fourier_terms(expression, t, specs)
    end
    @test_throws ArgumentError fourier_terms(cos(ω * t), t, [(ω, tan)])
    @test_throws ArgumentError fourier_terms(cos(ω * t), t, [(ω * t, cos)])
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

@testset "fraction-free linear solve" begin
    using LinearAlgebra
    using QuestBase: fraction_free_linear_solve

    @variables a b c x y z
    equations = [2x + y - z ~ a, x + 3y + z ~ b, 2x - y + 4z ~ c]
    solution = fraction_free_linear_solve(equations, [x, y, z])
    @eqtest solution[1] == (13a - 3b + 4c) / 31
    @eqtest solution[2] == (-2a + 10b - 3c) / 31
    @eqtest solution[3] == (-7a + 4b + 5c) / 31

    # A zero diagonal requires a row pivot.
    pivoted = fraction_free_linear_solve([y ~ a, x + y ~ b], [x, y])
    @eqtest pivoted[1] == b - a
    @eqtest pivoted[2] == a

    @test_throws LinearAlgebra.SingularException fraction_free_linear_solve(
        [x + y ~ a, 2x + 2y ~ b], [x, y]
    )

    # Compare arbitrary nonsingular integer systems against exact rational arithmetic.
    rng = Random.MersenneTwister(0x5eed)
    for dimension in 2:5
        matrix = rand(rng, -4:4, dimension, dimension)
        while iszero(det(matrix))
            matrix = rand(rng, -4:4, dimension, dimension)
        end
        rhs = rand(rng, -4:4, dimension)
        symbolic_variables = only(@variables q[1:dimension])
        random_equations = [
            sum(
                matrix[row, column] * symbolic_variables[column] for column in 1:dimension
            ) ~ rhs[row] for row in 1:dimension
        ]
        actual = fraction_free_linear_solve(random_equations, symbolic_variables)
        expected = Rational{Int}.(matrix) \ Rational{Int}.(rhs)
        for index in eachindex(actual)
            @eqtest actual[index] == expected[index]
        end
    end
end

@testset "fraction-free linear solve splits independent blocks" begin
    using LinearAlgebra
    using QuestBase: fraction_free_linear_solve, system_blocks

    @variables a b c x y z w p q

    # Two 2x2 blocks, interleaved so the decomposition cannot rely on contiguous indices.
    # Each solution carries only its own block's determinant: solved whole, every one of
    # them would come out over the product of both.
    equations = [p * x + q * z ~ a, x - z ~ b, p * y + w ~ c, y - w ~ 0]
    solution = fraction_free_linear_solve(equations, [x, y, z, w])
    @eqtest solution[1] == (a + b * q) / (p + q)
    @eqtest solution[3] == (a - b * p) / (p + q)
    @eqtest solution[2] == c / (1 + p)
    @eqtest solution[4] == c / (1 + p)

    # Transcendental coefficients do not change the decomposition, and the trig-free block
    # is not dragged through the rotation block's determinant.
    trig = fraction_free_linear_solve(
        [cos(p) * x - sin(p) * z ~ a, sin(p) * x + cos(p) * z ~ b, y ~ c], [x, y, z]
    )
    @eqtest trig[2] == c
    # Put the rotation on a rational point of the unit circle so the residual is exact.
    numeric = Dict(cos(p) => 3 // 5, sin(p) => 4 // 5, a => 1, b => 2)
    @eqtest substitute(cos(p) * trig[1] - sin(p) * trig[3], numeric) == 1
    @eqtest substitute(sin(p) * trig[1] + cos(p) * trig[3], numeric) == 2

    # Whole matrix nonzero: one component, so there is nothing to split and the solve runs
    # as before.
    @test system_blocks(Num[1 1; 1 1]) == [([1, 2], [1, 2])]

    # A structurally singular pattern has an unbalanced component, and is left to the full
    # solve to reject rather than being silently split.
    @test system_blocks(Num[1 1; 0 0]) == [([1, 2], [1, 2])]
    @test_throws LinearAlgebra.SingularException fraction_free_linear_solve(
        [x + y ~ a, 0 * x + 0 * y ~ b], [x, y]
    )

    # Random block-diagonal integer systems against exact rational arithmetic.
    rng = Random.MersenneTwister(0xb10c)
    for _ in 1:8
        sizes = rand(rng, 1:3, 3)
        dimension = sum(sizes)
        matrix = zeros(Int, dimension, dimension)
        offset = 0
        for size in sizes
            block = rand(rng, -4:4, size, size)
            while iszero(det(block))
                block = rand(rng, -4:4, size, size)
            end
            matrix[(offset + 1):(offset + size), (offset + 1):(offset + size)] = block
            offset += size
        end
        permutation = Random.randperm(rng, dimension)
        matrix = matrix[permutation, :]
        rhs = rand(rng, -4:4, dimension)
        symbolic_variables = only(@variables r[1:dimension])
        random_equations = [
            sum(
                matrix[row, column] * symbolic_variables[column] for column in 1:dimension
            ) ~ rhs[row] for row in 1:dimension
        ]
        actual = fraction_free_linear_solve(random_equations, symbolic_variables)
        expected = Rational{Int}.(matrix) \ Rational{Int}.(rhs)
        for index in eachindex(actual)
            @eqtest actual[index] == expected[index]
        end
    end
end

@testset "rearrange! reduces determinant denominators" begin
    # Fraction-free elimination produces one determinant denominator instead of nested
    # LU fractions. For a trigonometric ansatz that determinant contains
    # cos² + sin² = 1, which `rearrange!` must collapse before downstream averaging.
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
