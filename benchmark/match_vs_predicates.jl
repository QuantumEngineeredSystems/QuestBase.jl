"""
Benchmark: Moshi @match vs predicate chains (isadd/ismul/isdiv/ispow)

Tests all QuestBase functions that use predicate-based dispatch on BasicSymbolic variants.
Compares current implementation against @match-based alternatives.
"""

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
    _simplify_exp_products_mul,
    exp_to_trig,
    _exp_to_trig_node,
    _has_negative_coefficient,
    is_function,
    power_of,
    is_trig
using Symbolics: Symbolics, @variables, unwrap, expand, Num, Equation, Differential, substitute
using SymbolicUtils:
    SymbolicUtils,
    BasicSymbolic,
    Postwalk,
    isterm,
    ispow,
    isadd,
    isdiv,
    ismul,
    issym,
    add_with_div,
    is_literal_number,
    unwrap_const,
    sorted_arguments,
    arguments,
    operation

using Moshi.Match: @match
const BSImpl = SymbolicUtils.BasicSymbolicImpl
using SymbolicUtils: AddMulVariant

# ============================================================================
# @match-based alternatives
# ============================================================================

# --- 1. _apply_termwise (3 branches: isadd, ismul, isdiv) ---
function _apply_termwise_match(f, x::BasicSymbolic)
    @match x begin
        BSImpl.AddMul(; variant, args) => if variant == AddMulVariant.ADD
            sum(f(arg) for arg in args)
        else  # MUL
            prod(f(arg) for arg in args)
        end
        BSImpl.Div(; num, den) =>
            _apply_termwise_match(f, num) / _apply_termwise_match(f, den)
        _ => f(x)
    end
end

# --- 2. expand_exp_power (3 branches: isadd, ismul, ispow+isexp) ---
function expand_exp_power_match(expr::BasicSymbolic)
    @match expr begin
        BSImpl.AddMul(; variant, args) => if variant == AddMulVariant.ADD
            sum(expand_exp_power_match(arg) for arg in args)
        else  # MUL
            prod(expand_exp_power_match(arg) for arg in args)
        end
        BSImpl.Term(; f, args) && if f === (^) end =>
            let base = args[1]
                if isterm(base) && operation(base) === exp
                    exp(arguments(base)[1] * args[2])
                else
                    expr
                end
            end
        _ => expr
    end
end
expand_exp_power_match(expr) = expr

# --- 3. simplify_exp_products (4 branches: isadd, isdiv, ismul, ispow+isexp) ---
function simplify_exp_products_match(expr::BasicSymbolic)
    @match expr begin
        BSImpl.AddMul(; variant) => if variant == AddMulVariant.ADD
            _apply_termwise_match(simplify_exp_products_match, expr)
        else  # MUL
            _simplify_exp_products_mul(expr)
        end
        BSImpl.Div() =>
            _apply_termwise_match(simplify_exp_products_match, expr)
        BSImpl.Term(; f, args) && if f === (^) end =>
            let base = args[1]
                if isterm(base) && operation(base) === exp
                    exp(arguments(base)[1] * args[2])
                else
                    expr
                end
            end
        _ => expr
    end
end
simplify_exp_products_match(x) = x

# --- 4. get_independent (5 branches: isadd, ismul, isdiv, ispow, isterm/issym) ---
function get_independent_match(x::BasicSymbolic, t::Num)
    @match x begin
        BSImpl.AddMul(; variant, args) => if variant == AddMulVariant.ADD
            sum(get_independent_match(arg, t) for arg in args)
        else  # MUL
            prod(get_independent_match(arg, t) for arg in args)
        end
        BSImpl.Div(; num, den) =>
            !is_function(den, t) ? get_independent_match(num, t) / den : 0
        BSImpl.Term(; f, args) => if f === (^)
            !is_function(args[1], t) && !is_function(args[2], t) ? x : 0
        else
            !is_function(x, t) ? x : 0
        end
        BSImpl.Sym() => !is_function(x, t) ? x : 0
        _ => x
    end
end
get_independent_match(x, t::Num) = x

# --- 5. _get_all_terms (3 branches: isadd, ismul, isdiv) ---
function _get_all_terms_match(x::BasicSymbolic)
    @match x begin
        BSImpl.AddMul(; variant, args) => if variant == AddMulVariant.ADD
            vcat([_get_all_terms_match(arg) for arg in args]...)
        else  # MUL
            collect(args)
        end
        BSImpl.Div(; num, den) =>
            [_get_all_terms_match(num)..., _get_all_terms_match(den)...]
        _ => [x]
    end
end
_get_all_terms_match(x) = x

# --- 6. expand_fraction (2 branches: isadd/ismul, isdiv) ---
function expand_fraction_match(x::BasicSymbolic)
    @match x begin
        BSImpl.AddMul() =>
            _apply_termwise_match(expand_fraction_match, x)
        BSImpl.Div(; num, den) => begin
            num_expanded = SymbolicUtils.expand(num)
            if isadd(num_expanded)
                sum(expand_fraction_match(arg / den) for arg in arguments(num_expanded))
            else
                x
            end
        end
        _ => x
    end
end

# --- 7. simplify_complex (2 branches: isadd/ismul/isdiv vs leaf) ---
function simplify_complex_match(x::BasicSymbolic)
    @match x begin
        BSImpl.AddMul() => _apply_termwise_match(simplify_complex_match, x)
        BSImpl.Div() => _apply_termwise_match(simplify_complex_match, x)
        _ => begin
            v = unwrap_const(x)
            if v isa Complex && iszero(imag(v))
                return real(v)
            end
            x
        end
    end
end
simplify_complex_match(x::Complex) = isequal(x.im, 0) ? x.re : x.re + im * x.im
simplify_complex_match(x) = x

# --- 8. exp_to_trig (2 branches: isadd/isdiv/ismul vs leaf) ---
function exp_to_trig_match(x::BasicSymbolic)
    @match x begin
        BSImpl.AddMul() => _apply_termwise_match(exp_to_trig_match, x)
        BSImpl.Div() => _apply_termwise_match(exp_to_trig_match, x)
        _ => begin
            result = _exp_to_trig_node(x)
            isnothing(result) ? x : result
        end
    end
end

# --- 9. power_of (2 branches: ispow+issym, issym+issym) ---
function power_of_match(x::BasicSymbolic, y::BasicSymbolic)
    @match x begin
        BSImpl.Term(; f, args) && if f === (^) && issym(y) end =>
            isequal(args[1], y) ? unwrap_const(args[2]) : 0
        BSImpl.Sym() && if issym(y) end =>
            isequal(x, y) ? 1 : 0
        _ => 0
    end
end
power_of_match(x, y) = 0

# --- 10. is_trig (1 branch: ispow) ---
function is_trig_match(f::BasicSymbolic)
    @match f begin
        BSImpl.Term(; f=op, args) && if op === (^) end =>
            let base = args[1]
                isterm(base) && operation(base) ∈ [cos, sin]
            end
        BSImpl.Term(; f=op) =>
            op ∈ [cos, sin]
        _ => false
    end
end
is_trig_match(f) = false

# ============================================================================
# Setup test expressions
# ============================================================================
@variables t x(t) y(t) ω0 ω F k a b c f θ

# Various expression types for thorough testing
expr_add = unwrap(a^2 + b * c + a * b^2 + c^3 + a * b * c)  # Add
expr_mul = unwrap(a * b * c * (a + b) * (b + c))              # Mul
expr_div = unwrap((a^2 + b * c + c^3) / (a * b + c^2))       # Div
expr_pow = unwrap(a^3)                                         # Pow
expr_sym = unwrap(a)                                           # Sym

# Exponential expressions
expr_exp_pow = unwrap(exp(a)^3 + exp(b)^2 + a * exp(c)^4)
expr_exp_prod = unwrap(exp(3a) * exp(4b) + exp(a) * exp(c))

# Independent expressions
expr_indep_simple = unwrap(a^2 + b * c + a * b)
expr_indep_trig = unwrap(cos(f * t)^2 + a * sin(f * t) + b)

# Trig expressions
exp_sum = unwrap(exp(im * a) + exp(-im * a) + exp(im * (a + b)))

# Fraction expression
expr_frac = unwrap((a^2 + b * c + a * b^2 + c^3) / (a + b))

# ============================================================================
# Run benchmarks
# ============================================================================
println("=" ^ 80)
println("Moshi @match vs Predicate Chains Benchmark")
println("=" ^ 80)

function compare(name, f_pred, f_match, args...)
    print("\n### $name\n")
    t_pred = @benchmark $f_pred($(args)...)
    t_match = @benchmark $f_match($(args)...)
    med_pred = median(t_pred).time
    med_match = median(t_match).time
    ratio = med_pred / med_match
    status = ratio > 1.05 ? "@match FASTER" : ratio < 0.95 ? "predicates FASTER" : "~same"
    println("  predicates: $(round(med_pred, digits=1)) ns")
    println("  @match:     $(round(med_match, digits=1)) ns")
    println("  ratio:      $(round(ratio, digits=2))x  ($status)")
    return (name=name, pred_ns=med_pred, match_ns=med_match, ratio=ratio)
end

results = []

# 1. _apply_termwise
push!(results, compare("_apply_termwise (add)", _apply_termwise, _apply_termwise_match, simplify_complex, expr_add))
push!(results, compare("_apply_termwise (mul)", _apply_termwise, _apply_termwise_match, simplify_complex, expr_mul))
push!(results, compare("_apply_termwise (div)", _apply_termwise, _apply_termwise_match, simplify_complex, expr_div))

# 2. expand_exp_power
push!(results, compare("expand_exp_power", expand_exp_power, expand_exp_power_match, expr_exp_pow))

# 3. simplify_exp_products
push!(results, compare("simplify_exp_products", simplify_exp_products, simplify_exp_products_match, expr_exp_prod))

# 4. get_independent (most branches — expected biggest win)
push!(results, compare("get_independent (simple)", get_independent, get_independent_match, expr_indep_simple, t))
push!(results, compare("get_independent (trig)", get_independent, get_independent_match, expr_indep_trig, t))

# 5. _get_all_terms
push!(results, compare("_get_all_terms", QuestBase._get_all_terms, _get_all_terms_match, expr_add))

# 6. expand_fraction
push!(results, compare("expand_fraction", expand_fraction, expand_fraction_match, expr_frac))

# 7. simplify_complex
push!(results, compare("simplify_complex (add)", simplify_complex, simplify_complex_match, expr_add))

# 8. exp_to_trig
push!(results, compare("exp_to_trig", exp_to_trig, exp_to_trig_match, exp_sum))

# 9. power_of
push!(results, compare("power_of", power_of, power_of_match, expr_pow, unwrap(a)))

# 10. is_trig
push!(results, compare("is_trig", is_trig, is_trig_match, unwrap(cos(f * t))))

# ============================================================================
# Summary table
# ============================================================================
println("\n\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)
println()
println("| Function | Predicates (ns) | @match (ns) | Ratio | Winner |")
println("|---|---|---|---|---|")
for r in results
    status = r.ratio > 1.05 ? "@match" : r.ratio < 0.95 ? "predicates" : "~same"
    println("| $(r.name) | $(round(r.pred_ns, digits=1)) | $(round(r.match_ns, digits=1)) | $(round(r.ratio, digits=2))x | $status |")
end
