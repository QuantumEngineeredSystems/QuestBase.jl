using Test
using Symbolics
using QuestBase: QuestBase, @eqtest

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

@testset "trig_reduce linearizes products of trig functions" begin
    using QuestBase: trig_reduce
    @variables a b x y

    @eqtest trig_reduce(cos(x)^2) == 1//2 + cos(2x) * (1//2)
    @eqtest trig_reduce(sin(x)^2) == 1//2 - cos(2x) * (1//2)
    @eqtest trig_reduce(2 * sin(x) * cos(x)) == sin(2x)
    # two independent arguments combine into sum and difference angles
    product = trig_reduce(cos(x) * cos(y))
    difference = product - (cos(x + y) + cos(x - y)) / 2
    @test Symbolics.value(
        substitute(difference, Dict(x => 0.3, y => 2.0); fold=Val(true))
    ) ≈ 0 atol = 1e-14
    # coefficients and t-free factors ride along
    @eqtest trig_reduce(a * cos(x)^2 + a * sin(x)^2 + b) == a + b
    # a non-polynomial trig factor is left as it stands
    @eqtest trig_reduce(a / cos(x)) == a / cos(x)
end

@testset "trig_reduce agrees numerically with its input" begin
    using QuestBase: trig_reduce
    @variables a b x y
    numeric(e, point) = Symbolics.value(substitute(e, point; fold=Val(true)))
    point = Dict(a => 0.7, b => -1.3, x => 0.41, y => 1.9)
    for e in [
        (a * cos(x) + b * sin(2x))^3,
        cos(x + y)^2 * sin(x - 2y),
        (a + cos(x)^2) * (b + sin(y)^3),
        cos(x)^4 + 2 * cos(x)^2 * sin(x)^2 + sin(x)^4,
    ]
        @test numeric(trig_reduce(e), point) ≈ numeric(e, point) atol = 1e-12
    end
end
