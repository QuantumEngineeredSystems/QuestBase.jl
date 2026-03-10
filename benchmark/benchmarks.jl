using BenchmarkTools
using QuestBase
using QuestBase:
    expand_all,
    expand_fraction,
    _apply_termwise,
    simplify_complex,
    get_independent,
    get_all_terms,
    expand_exp_power,
    simplify_exp_products,
    trig_to_exp,
    exp_to_trig,
    trig_reduce,
    fourier_cos_term,
    fourier_sin_term,
    drop_powers,
    max_power,
    power_of,
    substitute_all,
    is_harmonic,
    DifferentialEquation,
    HarmonicVariable,
    d
using Symbolics: Symbolics, @variables, unwrap, expand, Num, Equation
using SymbolicUtils: BasicSymbolic

const SUITE = BenchmarkGroup()

# --- Setup variables ---
@variables t x(t) y(t) ω0 ω F k a b c f θ

# --- Symbolic utilities ---
SUITE["symbolic_utils"] = BenchmarkGroup()

expr_add = a + b + c
expr_mul = a * b * c
expr_div = a / b
expr_nested = (a + b * c)^2 + (a - c) / b
expr_complex = a + b * cos(f * t) + c * sin(f * t)

SUITE["symbolic_utils"]["expand_all_simple"] = @benchmarkable expand_all($expr_add)
SUITE["symbolic_utils"]["expand_all_nested"] = @benchmarkable expand_all($expr_nested)
SUITE["symbolic_utils"]["expand_all_trig"] = @benchmarkable expand_all($expr_complex)

SUITE["symbolic_utils"]["expand_fraction"] = @benchmarkable expand_fraction(
    $(unwrap((a + b) / c))
)
SUITE["symbolic_utils"]["_apply_termwise_add"] = @benchmarkable _apply_termwise(
    identity, $(unwrap(expr_add))
)
SUITE["symbolic_utils"]["_apply_termwise_mul"] = @benchmarkable _apply_termwise(
    identity, $(unwrap(expr_mul))
)
SUITE["symbolic_utils"]["_apply_termwise_div"] = @benchmarkable _apply_termwise(
    identity, $(unwrap(expr_div))
)

SUITE["symbolic_utils"]["simplify_complex_add"] = @benchmarkable simplify_complex(
    $(unwrap(Complex{Num}(a + b, 0 * a)))
)
SUITE["symbolic_utils"]["simplify_complex_real"] = @benchmarkable simplify_complex(
    $(Complex{Num}(a, 0 * a))
)

SUITE["symbolic_utils"]["get_independent_simple"] = @benchmarkable get_independent(
    $a + $b + $c, $t
)
SUITE["symbolic_utils"]["get_independent_trig"] = @benchmarkable get_independent(
    $(cos(f * t)^2 + a + b), $t
)

SUITE["symbolic_utils"]["get_all_terms_add"] = @benchmarkable get_all_terms($expr_add)
SUITE["symbolic_utils"]["get_all_terms_nested"] = @benchmarkable get_all_terms($expr_nested)

# --- Exponentials ---
SUITE["exponentials"] = BenchmarkGroup()

SUITE["exponentials"]["expand_exp_power"] = @benchmarkable expand_exp_power(
    $(unwrap(exp(a)^3))
)
SUITE["exponentials"]["expand_exp_power_nested"] = @benchmarkable expand_exp_power(
    $(unwrap(exp(a)^3 + exp(b)^2))
)

SUITE["exponentials"]["simplify_exp_products"] = @benchmarkable simplify_exp_products(
    $(unwrap(exp(a) * exp(b)))
)
SUITE["exponentials"]["simplify_exp_products_multi"] = @benchmarkable simplify_exp_products(
    $(unwrap(exp(3a) * exp(4b)))
)

# --- Fourier ---
SUITE["fourier"] = BenchmarkGroup()

trig_expr_simple = cos(f * t)
trig_expr_sq = cos(f * t)^2
trig_expr_product =
    (a * sin(f * t) + b * cos(f * t)) *
    (a * sin(2 * f * t) + b * cos(2 * f * t)) *
    (a * sin(3 * f * t) + b * cos(3 * f * t))
trig_expr_hard = (a + b * cos(f * t + θ)^2)^3 * sin(f * t)

SUITE["fourier"]["trig_to_exp_simple"] = @benchmarkable trig_to_exp($trig_expr_simple)
SUITE["fourier"]["trig_to_exp_sq"] = @benchmarkable trig_to_exp($trig_expr_sq)

SUITE["fourier"]["exp_to_trig"] = @benchmarkable exp_to_trig($(unwrap(exp(im * a))))

SUITE["fourier"]["trig_reduce_simple"] = @benchmarkable trig_reduce($trig_expr_sq)
SUITE["fourier"]["trig_reduce_product"] = @benchmarkable trig_reduce($trig_expr_product)

SUITE["fourier"]["fourier_cos_term_simple"] = @benchmarkable fourier_cos_term(
    $trig_expr_simple, $f, $t
)
SUITE["fourier"]["fourier_cos_term_sq"] = @benchmarkable fourier_cos_term(
    $(cos(f * t)^2), $(2f), $t
)
SUITE["fourier"]["fourier_sin_term_simple"] = @benchmarkable fourier_sin_term(
    $(sin(f * t)), $f, $t
)
SUITE["fourier"]["fourier_cos_term_product"] = @benchmarkable fourier_cos_term(
    $trig_expr_product, $(2f), $t
)
SUITE["fourier"]["fourier_cos_term_hard"] = @benchmarkable fourier_cos_term(
    $trig_expr_hard, $f, $t
)

# --- Powers ---
SUITE["powers"] = BenchmarkGroup()

SUITE["powers"]["max_power_simple"] = @benchmarkable max_power($(a^2 + b), $a)
SUITE["powers"]["max_power_nested"] = @benchmarkable max_power($(a * ((a + b)^4)^2 + a), $a)
SUITE["powers"]["power_of_pow"] = @benchmarkable power_of($(unwrap(a^3)), $(unwrap(a)))
SUITE["powers"]["drop_powers_simple"] = @benchmarkable drop_powers($(a^2 + b), $a, 1)
SUITE["powers"]["drop_powers_multi"] = @benchmarkable drop_powers($((a + b)^2), [$a, $b], 1)
SUITE["powers"]["drop_powers_high"] = @benchmarkable drop_powers(
    $((a + b)^3 + (a + b)^5), [$a, $b], 4
)

# --- Substitution ---
SUITE["substitution"] = BenchmarkGroup()

@variables d_var e_var f_var g_var h_var
rules_small = Dict(a => b)
rules_medium = Dict(a => b, c => d_var, e_var => f_var)
expr_sub = a + b + c + d_var + e_var + f_var

SUITE["substitution"]["substitute_all_1rule"] = @benchmarkable substitute_all(
    $(a + b), $(a => b)
)
SUITE["substitution"]["substitute_all_3rules"] = @benchmarkable substitute_all(
    $expr_sub, $rules_medium
)
SUITE["substitution"]["substitute_all_vector"] = @benchmarkable substitute_all(
    [$a, $c, $e_var], $rules_medium
)

# --- Construction ---
SUITE["construction"] = BenchmarkGroup()

SUITE["construction"]["DifferentialEquation_single"] = @benchmarkable DifferentialEquation(
    d($x, $t, 2) + $ω0^2 * $x ~ $F * cos($ω * $t), $x
)
SUITE["construction"]["DifferentialEquation_coupled"] = @benchmarkable DifferentialEquation(
    [d($x, $t, 2) + $ω0^2 * $x - $k * $y, d($y, $t, 2) + $ω0^2 * $y - $k * $x] .~
    [$F * cos($ω * $t), 0],
    [$x, $y],
)

# --- Harmonic checks ---
SUITE["harmonic"] = BenchmarkGroup()

SUITE["harmonic"]["is_harmonic_true"] = @benchmarkable is_harmonic($(cos(f * t)), $t)
SUITE["harmonic"]["is_harmonic_false"] = @benchmarkable is_harmonic($(cos(f * t^2 + a)), $t)
