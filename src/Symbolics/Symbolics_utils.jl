"The vartype every expression here is built on, for callers with no `BasicSymbolic` to read it off."
const NUM_VARTYPE = vartype(fieldtype(Num, 1))

"The concrete expression node every real-valued `Num` wraps."
const SymbolicNode = BasicSymbolic{NUM_VARTYPE}
const SymbolicArguments = SymbolicUtils.ROArgsT{NUM_VARTYPE}
const SYMBOLIC_ONE = Const{NUM_VARTYPE}(1)
const SYMBOLIC_ZERO = Const{NUM_VARTYPE}(0)
const SYMBOLIC_TWO = Const{NUM_VARTYPE}(2)

"`x` as `base^exponent`, with exponent one for anything that is not a power."
function base_exponent(x::SymbolicNode)
    return @match x begin
        BSImpl.Term(; f, args) && if f === (^)
        end => (args[1], args[2])
        _ => (x, SYMBOLIC_ONE)
    end
end

"`x` as `numerator / denominator`, with denominator one for anything that is not a fraction."
function numerator_denominator(x::SymbolicNode)
    return @match x begin
        BSImpl.Div(; num, den) => (num, den)
        _ => (x, SYMBOLIC_ONE)
    end
end

"Return true if `x` is a call to `cos` or `sin`."
function is_trig_call(x::SymbolicNode)
    return @match x begin
        BSImpl.Term(; f) => f === cos || f === sin
        _ => false
    end
end

_apply_termwise(f, x::Num) = Num(_apply_termwise(f, unwrap(x)))

function expand_fraction(x::SymbolicNode)
    return @match x begin
        BSImpl.AddMul(;) => _apply_termwise(expand_fraction, x)
        BSImpl.Div(; num, den) => begin
            expanded = SymbolicUtils.expand(num)
            isadd(expanded) || return x
            add_worker(
                NUM_VARTYPE,
                SymbolicNode[expand_fraction(arg / den) for arg in arguments(expanded)],
            )
        end
        _ => x
    end
end
expand_fraction(x::Num) = Num(expand_fraction(unwrap(x)))

"""
    with_number(f, value)

`f` applied to the number `value` with its concrete type. SymbolicUtils stores constants and
the coefficients of sums and products untyped; splitting them over the literal types here
keeps every reader inferable.
"""
function with_number(f::F, @nospecialize(value)) where {F}
    return @match value begin
        number::Int => f(number)
        number::Float64 => f(number)
        number::Rational{Int} => f(number)
        number::BigInt => f(number)
        number::Rational{BigInt} => f(number)
        number::BigFloat => f(number)
        number::Float32 => f(number)
        number::Bool => f(number)
        number::Complex{Int} => f(number)
        number::ComplexF64 => f(number)
        number::Complex{Bool} => f(number)
        number::Complex{Rational{Int}} => f(number)
        number::Irrational{:π} => f(number)
        number::Irrational{:ℯ} => f(number)
        number::Irrational{:γ} => f(number)
        number::Irrational{:φ} => f(number)
        number::Irrational{:catalan} => f(number)
        _ => throw(ArgumentError("Unsupported constant $value of type $(typeof(value))."))
    end
end

"`f` applied to the payload of the constant `x` with its concrete type."
with_const(f::F, x::SymbolicNode) where {F} = with_number(f, unwrap_const(x))

"`f` applied to the payload of the constant `x`, which must be an exact rational."
function with_exact(f::F, x::SymbolicNode) where {F}
    return @match unwrap_const(x) begin
        value::Int => f(value)
        value::Rational{Int} => f(value)
        value::BigInt => f(value)
        value::Rational{BigInt} => f(value)
        value::Bool => f(value)
        _ => throw(ArgumentError("$x is not an exact rational constant."))
    end
end

"""
Apply a function f on every member of a sum or a product.

Rebuilt with `add_worker`/`mul_worker`; folding with `sum` is `O(n²)`.
"""
function _apply_termwise(f::F, x::SymbolicNode) where {F}
    return @match x begin
        BSImpl.AddMul(; variant) => begin
            terms = SymbolicNode[f(arg) for arg in arguments(x)]
            if variant === AddMulVariant.ADD
                add_worker(NUM_VARTYPE, terms)
            else
                mul_worker(NUM_VARTYPE, terms)
            end
        end
        BSImpl.Div(; num, den) => _apply_termwise(f, num) / _apply_termwise(f, den)
        _ => f(x)
    end
end

simplify_complex(x::Complex) = isequal(x.im, 0) ? x.re : x.re + im * x.im
simplify_complex(x) = x
function simplify_complex(x::SymbolicNode)
    return @match x begin
        BSImpl.AddMul(;) || BSImpl.Div(;) => _apply_termwise(simplify_complex, x)
        BSImpl.Const(;) =>
            with_const(has_zero_imaginary, x) ? with_const(real_constant, x) : x
        _ => x
    end
end

"Return true if the complex constant `value` is real."
has_zero_imaginary(value::Complex) = iszero(imag(value))
has_zero_imaginary(value::Real) = false

"The real part of `value` as a symbolic constant."
real_constant(value::Number) = Const{NUM_VARTYPE}(real(value))

const Subtype = Union{Num,Equation,BasicSymbolic}

function _deep_substitute(e::BasicSymbolic, unwrap_rules::Dict)
    return SymbolicUtils.substitute(e, unwrap_rules; filterer=_ -> true)
end
_deep_substitute(x::Num, ur::Dict) = wrap(_deep_substitute(unwrap(x), ur))
function _deep_substitute(eq::Equation, ur::Dict)
    return Equation(
        _deep_substitute(unwrap(eq.lhs), ur), _deep_substitute(unwrap(eq.rhs), ur)
    )
end
_deep_substitute(x, ::Dict) = x

"""
$(TYPEDSIGNATURES)

Perform substitutions in `rules` on `x`.
`include_derivatives=true` also substitutes inside derivative arguments using
`Symbolics.substitute_in_deriv`.
"""
function substitute_all(x::Subtype, rules::Dict; include_derivatives=true)
    unwrap_rules = Dict(unwrap(k) => unwrap(v) for (k, v) in rules)
    result = _deep_substitute(x, unwrap_rules)
    if include_derivatives
        result = Symbolics.substitute_in_deriv(result, rules)
    end
    return result
end

const Collections = Union{Dict,Pair,Vector,OrderedDict}
substitute_all(v::AbstractArray, rules) = [substitute_all(x, rules) for x in v]
substitute_all(x::Subtype, rules::Collections) = substitute_all(x, Dict(rules))
function substitute_all(x::Complex{Num}, rules::Collections)
    return substitute_all(x.re, rules) + im * substitute_all(x.im, rules)
end

get_independent(x::Num, t::Num) = Num(get_independent(unwrap(x), t))
function get_independent(x::Complex{Num}, t::Num)
    return get_independent(x.re, t) + im * get_independent(x.im, t)
end
get_independent(v::Vector{Num}, t::Num) = [get_independent(el, t) for el in v]
get_independent(x, t::Num) = x

function get_independent(x::SymbolicNode, t::Num)
    return @match x begin
        BSImpl.AddMul(;) => _apply_termwise(arg -> get_independent(arg, t), x)
        BSImpl.Div(; num, den) =>
            is_function(den, t) ? SYMBOLIC_ZERO : get_independent(num, t) / den
        BSImpl.Term(;) || BSImpl.Sym(;) => is_function(x, t) ? SYMBOLIC_ZERO : x
        _ => x
    end
end

"Return all the terms contained in `x`"
get_all_terms(x::Num) = Num.(get_all_terms(unwrap(x)))
get_all_terms(x::SymbolicNode) = unique(append_terms!(SymbolicNode[], Symbolics.expand(x)))
function get_all_terms(x::Equation)
    return unique(cat(get_all_terms(Num(x.lhs)), get_all_terms(Num(x.rhs)); dims=1))
end

"Push the summands of `x` onto `terms`, split into their factors and fraction parts."
function append_terms!(terms::Vector{SymbolicNode}, x::SymbolicNode)
    @match x begin
        BSImpl.AddMul(; variant) => if variant === AddMulVariant.ADD
            for arg in arguments(x)
                append_terms!(terms, arg)
            end
        else
            append!(terms, arguments(x))
        end
        BSImpl.Div(; num, den) => begin
            append_terms!(terms, num)
            append_terms!(terms, den)
        end
        _ => push!(terms, x)
    end
    return terms
end

function is_harmonic(x::Num, t::Num)
    all_terms = get_all_terms(Num(time_dependent_part(unwrap(x), t)))
    t_terms = setdiff(all_terms, get_independent(all_terms, t))
    return all(term -> is_harmonic_factor(unwrap(term), t), t_terms)
end

"""
Return true if `term` is a positive integer power of a cosine or sine whose argument is
linear in `t`. A power such as `cos(ωt)^2` is a finite Fourier series too, and the
projection expands it exactly.
"""
function is_harmonic_factor(term::SymbolicNode, t::Num)
    base, exponent = base_exponent(term)
    is_trig_call(base) || return false
    (isconst(exponent) && with_const(is_expanded_power, exponent)) || return false
    return isone(max_power(Num(first(arguments(base))), t))
end

"""
The summands of `x` that are functions of `t`. Expanding a summand free of `t` yields only
terms free of `t`, and those can neither merge with nor cancel a term holding `t`, so the
terms of `x` that depend on `t` are those of this part. It is far smaller to expand.
"""
function time_dependent_part(x::SymbolicNode, t::Num)
    return @match x begin
        BSImpl.AddMul(; dict, variant) && if variant === AddMulVariant.ADD
        end => SymbolicUtils.Add{NUM_VARTYPE}(
            0,
            Dict{SymbolicNode,Number}(
                summand => weight for (summand, weight) in dict if is_function(summand, t)
            );
            type=Real,
        )
        _ => x
    end
end

is_harmonic(x::Equation, t::Num) = is_harmonic(x.lhs, t) && is_harmonic(x.rhs, t)
is_harmonic(x, t) = is_harmonic(Num(x), Num(t))

"Return true if `f` is a function of `var`."
is_function(f, var) = unwrap(var) in get_variables(f)

"""
$(TYPEDSIGNATURES)

Return true if `needle` occurs as a subexpression of `haystack`.

Structural, not `occursin(string(needle), string(haystack))`: pretty-printing costs a full
render and normalises, reorders and elides.
"""
occurs_in(needle, haystack) = _occurs_in(unwrap(needle), unwrap(haystack))

function _occurs_in(target, expr)
    isequal(expr, target) && return true
    expr isa BasicSymbolic && SymbolicUtils.iscall(expr) || return false
    return any(arg -> _occurs_in(target, arg), arguments(expr))
end

"""
$(TYPEDSIGNATURES)

Return the subexpressions of `x` in which a variable of `vars` (or a derivative of one)
appears inside an operation other than an addition, a multiplication or a power with a
non-negative integer exponent. An empty result means `x` is a polynomial in `vars`.

Averaging replaces each variable by a truncated Fourier series and projects the result onto
the harmonics of the ansatz. Only polynomials map a finite Fourier series onto a finite
Fourier series, so terms such as `sin(x)`, `exp(x)` or `1/x` cannot be averaged and are
returned here.
"""
function nonpolynomial_terms(x, vars)
    unwrapped = vars isa AbstractVector ? unwrap.(vars) : [unwrap(vars)]
    return unique(collect_nonpolynomial!(BasicSymbolic[], unwrap(x), unwrapped))
end
function nonpolynomial_terms(eq::Equation, vars)
    return unique(
        cat(nonpolynomial_terms(eq.lhs, vars), nonpolynomial_terms(eq.rhs, vars); dims=1)
    )
end

"""
$(TYPEDSIGNATURES)

Return true if `x` is a polynomial in `vars` and their derivatives, i.e. if it has no
[`nonpolynomial_terms`](@ref).
"""
is_polynomial(x, vars) = isempty(nonpolynomial_terms(x, vars))

"Return true if `x` is one of `vars` or contains one of them as a subexpression."
contains_variable(x, vars) = any(v -> isequal(x, v), vars)
function contains_variable(x::BasicSymbolic, vars)
    any(v -> isequal(x, v), vars) && return true
    SymbolicUtils.iscall(x) || return false
    return any(arg -> contains_variable(arg, vars), arguments(x))
end

"Walk `x` and push each subexpression which is not polynomial in `vars` onto `terms`."
function collect_nonpolynomial!(terms, x, vars)
    contains_variable(x, vars) || return terms
    any(v -> isequal(x, v), vars) && return terms

    @match x begin
        BSImpl.AddMul(;) => for arg in arguments(x)
            collect_nonpolynomial!(terms, arg, vars)
        end
        BSImpl.Div(; num, den) => begin
            collect_nonpolynomial!(terms, num, vars)
            contains_variable(den, vars) && push!(terms, x)
        end
        BSImpl.Term(; f, args) => if f === (^)
            base, exponent = args
            if isconst(exponent) && with_const(is_natural_exponent, exponent)
                collect_nonpolynomial!(terms, base, vars)
            else
                push!(terms, x)
            end
        elseif f isa Differential
            collect_nonpolynomial!(terms, first(args), vars)
        else
            push!(terms, x)
        end
        _ => push!(terms, x)
    end
    return terms
end

"Return true if `power` is a non-negative integer, so that raising to it keeps a polynomial."
is_natural_exponent(power::Real) = isinteger(power) && power >= 0
is_natural_exponent(power::Complex) = false

"""
Counts the number of derivatives of a symbolic variable.
"""
function count_derivatives(x::SymbolicNode)
    return @match x begin
        BSImpl.Term(; f) => f isa Differential ? f.order : 0
        BSImpl.Sym(;) => 0
        _ => throw(ArgumentError("The input is not a single term or symbol"))
    end
end
count_derivatives(x::Num) = count_derivatives(unwrap(x))
