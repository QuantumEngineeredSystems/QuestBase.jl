using Test
using Symbolics
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
