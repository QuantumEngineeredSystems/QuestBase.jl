using Test
using Symbolics
using QuestBase: QuestBase, @eqtest

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

@testset "fourier_series returns every harmonic" begin
    using QuestBase: fourier_series, fourier_terms
    @variables ω t a b θ
    x = (a * cos(ω * t) + b * sin(ω * t))^3 + 2 + a * cos(ω * t + θ)^2
    series = fourier_series(x, t)
    @eqtest getfield.(series, :frequency) == [0, ω, 2ω, 3ω]
    for component in series
        specs = [(component.frequency, cos), (component.frequency, sin)]
        cosine, sine = fourier_terms(x, t, specs)
        @eqtest component.cosine == cosine
        isequal(component.frequency, 0) || @eqtest component.sine == sine
    end
    @test isempty(fourier_series(Num(0), t))
end
