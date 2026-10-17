using Test
using Symbolics
using QuestBase: QuestBase, @eqtest

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
