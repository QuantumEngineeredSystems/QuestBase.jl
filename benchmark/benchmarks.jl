using BenchmarkTools
using QuestBase
using QuestBase:
    expand_fraction,
    _apply_termwise,
    simplify_complex,
    get_independent,
    get_all_terms,
    trig_reduce,
    fourier_cos_term,
    fourier_sin_term,
    fourier_series,
    fraction_free_linear_solve,
    drop_powers,
    max_power,
    power_of,
    substitute_all,
    is_harmonic,
    is_function,
    count_derivatives,
    DifferentialEquation,
    HarmonicVariable,
    HarmonicEquation,
    add_harmonic!,
    rearrange_standard!,
    d,
    var_name
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

SUITE["symbolic_utils"]["simplify_complex_add"] = @benchmarkable simplify_complex(
    $(unwrap(Complex{Num}(a^2 + b * c, 0 * a)))
)
SUITE["symbolic_utils"]["simplify_complex_real"] = @benchmarkable simplify_complex(
    $(Complex{Num}(a^2 + b * c, 0 * a))
)

SUITE["symbolic_utils"]["get_independent_simple"] = @benchmarkable get_independent(
    $(a^2 + b * c + a * b), $t
)
SUITE["symbolic_utils"]["get_independent_trig"] = @benchmarkable get_independent(
    $(cos(f * t)^2 + a * sin(f * t) + b), $t
)

SUITE["symbolic_utils"]["get_all_terms_nested"] = @benchmarkable get_all_terms($expr_nested)

SUITE["symbolic_utils"]["is_function_true"] = @benchmarkable is_function(
    $(cos(f * t)^2 + a), $t
)
SUITE["symbolic_utils"]["is_function_false"] = @benchmarkable is_function(
    $(a^2 + b * c), $t
)

SUITE["symbolic_utils"]["count_derivatives_0"] = @benchmarkable count_derivatives($x)
SUITE["symbolic_utils"]["count_derivatives_2"] = @benchmarkable count_derivatives(
    $(d(d(x, t), t))
)

SUITE["symbolic_utils"]["var_name"] = @benchmarkable var_name($x)

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

SUITE["fourier"]["trig_reduce_hard"] = @benchmarkable trig_reduce($trig_expr_hard)
SUITE["fourier"]["fourier_series_product"] = @benchmarkable fourier_series(
    $trig_expr_product, $t
)
SUITE["fourier"]["fourier_series_hard"] = @benchmarkable fourier_series($trig_expr_hard, $t)

# ==========================================================================
# Linear solve
# ==========================================================================
SUITE["linear_solve"] = BenchmarkGroup()

@variables u1 v1 u2 v2 γ η
dense_targets = [u1, v1, u2]
dense_system = [
    (a + b * c) * u1 + (b - γ) * v1 + c * u2 ~ F,
    (c * γ - a) * u1 + (a + b^2) * v1 + η * u2 ~ 0,
    (η * a) * u1 + (γ + c) * v1 + (a * b - η) * u2 ~ k,
]
SUITE["linear_solve"]["dense_3x3"] = @benchmarkable fraction_free_linear_solve(
    $dense_system, $dense_targets
)

rotation = [
    cos(ω * t) sin(ω * t) 0 0
    -sin(ω * t) cos(ω * t) 0 0
    0 0 cos(3ω * t) sin(3ω * t)
    0 0 -sin(3ω * t) cos(3ω * t)
]
mass_targets = [u1, v1, u2, v2]
mass_system = rotation * mass_targets .~ [a * u1 + γ * v1, b * v1 - u1, c * u2, η * v2 + F]
SUITE["linear_solve"]["rotating_mass_matrix"] = @benchmarkable fraction_free_linear_solve(
    $mass_system, $mass_targets
)

# ==========================================================================
# Powers
# ==========================================================================
SUITE["powers"] = BenchmarkGroup()

SUITE["powers"]["max_power_simple"] = @benchmarkable max_power($(a^2 + b), $a)
SUITE["powers"]["max_power_nested"] = @benchmarkable max_power($(a * ((a + b)^4)^2 + a), $a)
SUITE["powers"]["power_of_pow"] = @benchmarkable power_of($(unwrap(a^3)), $(unwrap(a)))
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
