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
    substitute_all,
    is_harmonic,
    count_derivatives,
    add_div,
    DifferentialEquation,
    HarmonicVariable,
    HarmonicEquation,
    add_harmonic!,
    rearrange_standard!,
    d
using Symbolics: Symbolics, @variables, unwrap, expand, Num, Equation, Differential
using SymbolicUtils: BasicSymbolic

const SUITE = BenchmarkGroup()

# --- Setup variables ---
@variables t x(t) y(t) ω0 ω F k a b c f θ
D = Differential(t)

# ==========================================================================
# Symbolic utilities
# ==========================================================================
SUITE["symbolic_utils"] = BenchmarkGroup()

expr_nested = (a + b * c)^2 + (a - c) / b
expr_complex = a + b * cos(f * t) + c * sin(f * t)

SUITE["symbolic_utils"]["expand_all_nested"] = @benchmarkable expand_all($expr_nested)
SUITE["symbolic_utils"]["expand_all_trig"] = @benchmarkable expand_all($expr_complex)

SUITE["symbolic_utils"]["expand_fraction"] = @benchmarkable expand_fraction(
    $(unwrap((a^2 + b * c + a * b^2 + c^3) / (a + b)))
)

SUITE["symbolic_utils"]["_apply_termwise_add"] = @benchmarkable _apply_termwise(
    simplify_complex, $(unwrap(a^2 + b * c + a * b^2 + c^3 + a * b * c))
)
SUITE["symbolic_utils"]["_apply_termwise_mul"] = @benchmarkable _apply_termwise(
    simplify_complex, $(unwrap(a * b * c * (a + b) * (b + c)))
)
SUITE["symbolic_utils"]["_apply_termwise_div"] = @benchmarkable _apply_termwise(
    simplify_complex, $(unwrap((a^2 + b * c + c^3) / (a * b + c^2)))
)

SUITE["symbolic_utils"]["simplify_complex"] = @benchmarkable simplify_complex(
    $(Complex{Num}(a^2 + b * c, 0 * a))
)

SUITE["symbolic_utils"]["get_independent_simple"] = @benchmarkable get_independent(
    $(a^2 + b * c + a * b), $t
)
SUITE["symbolic_utils"]["get_independent_trig"] = @benchmarkable get_independent(
    $(cos(f * t)^2 + a * sin(f * t) + b), $t
)

SUITE["symbolic_utils"]["get_all_terms_nested"] = @benchmarkable get_all_terms($expr_nested)

SUITE["symbolic_utils"]["count_derivatives"] = @benchmarkable count_derivatives(
    $(d(d(x, t), t))
)

# ==========================================================================
# Exponentials
# ==========================================================================
SUITE["exponentials"] = BenchmarkGroup()

SUITE["exponentials"]["expand_exp_power"] = @benchmarkable expand_exp_power(
    $(unwrap(exp(a)^3 + exp(b)^2 + a * exp(c)^4))
)
SUITE["exponentials"]["simplify_exp_products"] = @benchmarkable simplify_exp_products(
    $(unwrap(exp(3a) * exp(4b) + exp(a) * exp(c)))
)

# ==========================================================================
# Fourier
# ==========================================================================
SUITE["fourier"] = BenchmarkGroup()

trig_expr_sq = cos(f * t)^2
trig_expr_product =
    (a * sin(f * t) + b * cos(f * t)) *
    (a * sin(2 * f * t) + b * cos(2 * f * t)) *
    (a * sin(3 * f * t) + b * cos(3 * f * t))
trig_expr_hard = (a + b * cos(f * t + θ)^2)^3 * sin(f * t)

SUITE["fourier"]["trig_to_exp_sq"] = @benchmarkable trig_to_exp($trig_expr_sq)
SUITE["fourier"]["trig_to_exp_product"] = @benchmarkable trig_to_exp($trig_expr_product)

exp_sum = trig_to_exp(a * cos(f * t) + b * sin(2f * t))
SUITE["fourier"]["exp_to_trig"] = @benchmarkable exp_to_trig($exp_sum)

SUITE["fourier"]["add_div"] = @benchmarkable add_div($(a / (b + c) + b / (a + c)))

SUITE["fourier"]["trig_reduce_simple"] = @benchmarkable trig_reduce($trig_expr_sq)
SUITE["fourier"]["trig_reduce_product"] = @benchmarkable trig_reduce($trig_expr_product)

SUITE["fourier"]["fourier_cos_term_simple"] = @benchmarkable fourier_cos_term(
    $(cos(f * t)), $f, $t
)
SUITE["fourier"]["fourier_cos_term_sq"] = @benchmarkable fourier_cos_term(
    $trig_expr_sq, $(2f), $t
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

# ==========================================================================
# Powers
# ==========================================================================
SUITE["powers"] = BenchmarkGroup()

SUITE["powers"]["max_power_simple"] = @benchmarkable max_power($(a^2 + b), $a)
SUITE["powers"]["max_power_nested"] = @benchmarkable max_power($(a * ((a + b)^4)^2 + a), $a)
SUITE["powers"]["drop_powers_single_var"] = @benchmarkable drop_powers(
    $((a + b)^3 + a^2 * b + a * b^2), $a, 2
)
SUITE["powers"]["drop_powers_multi"] = @benchmarkable drop_powers($((a + b)^2), [$a, $b], 1)
SUITE["powers"]["drop_powers_high"] = @benchmarkable drop_powers(
    $((a + b)^3 + (a + b)^5), [$a, $b], 4
)

# ==========================================================================
# Substitution
# ==========================================================================
SUITE["substitution"] = BenchmarkGroup()

@variables d_var e_var f_var g_var h_var
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
SUITE["substitution"]["substitute_all_no_deriv"] = @benchmarkable substitute_all(
    $(a^2 + b * c), $(Dict(a => b)); include_derivatives=false
)
SUITE["substitution"]["substitute_all_with_deriv"] = @benchmarkable substitute_all(
    $(D(x) + a * x), $(Dict(x => 2x, a => b))
)

# ==========================================================================
# Construction and rearrangement
# ==========================================================================
SUITE["construction"] = BenchmarkGroup()

SUITE["construction"]["DifferentialEquation_single"] = @benchmarkable DifferentialEquation(
    d($x, $t, 2) + $ω0^2 * $x ~ $F * cos($ω * $t), $x
)
SUITE["construction"]["DifferentialEquation_coupled"] = @benchmarkable DifferentialEquation(
    [d($x, $t, 2) + $ω0^2 * $x - $k * $y, d($y, $t, 2) + $ω0^2 * $y - $k * $x] .~
    [$F * cos($ω * $t), 0],
    [$x, $y],
)

# Rearrangement benchmark
diff_eq_for_rearrange = DifferentialEquation(
    [d(x, t, 2) + ω0^2 * x - k * y ~ F * cos(ω * t), d(y, t, 2) + ω0^2 * y - k * x ~ 0],
    [x, y],
)
SUITE["construction"]["rearrange_standard"] = @benchmarkable rearrange_standard!(eom) setup = (
    eom = deepcopy($diff_eq_for_rearrange)
)

# ==========================================================================
# Harmonic checks
# ==========================================================================
SUITE["harmonic"] = BenchmarkGroup()

SUITE["harmonic"]["is_harmonic_true"] = @benchmarkable is_harmonic($(cos(f * t)), $t)
SUITE["harmonic"]["is_harmonic_false"] = @benchmarkable is_harmonic($(cos(f * t^2 + a)), $t)
SUITE["harmonic"]["is_harmonic_complex"] = @benchmarkable is_harmonic(
    $(a * cos(f * t) + b * sin(2f * t) + c), $t
)
