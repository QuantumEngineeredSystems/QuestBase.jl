
expand_all(x::Num) = Num(expand_all(unwrap(x)))
_apply_termwise(f, x::Num) = wrap(_apply_termwise(f, unwrap(x)))

"Expands using SymbolicUtils.expand and expand_exp_power (changes exp(x)^n to exp(x*n)"
function expand_all(x)
    result = Postwalk(expand_exp_power)(SymbolicUtils.expand(x))
    return isnothing(result) ? x : result
end
expand_all(x::Complex{Num}) = expand_all(x.re) + im * expand_all(x.im)

function expand_fraction(x::BasicSymbolic)
    if isadd(x) || ismul(x)
        _apply_termwise(expand_fraction, x)
    elseif isdiv(x)
        args = arguments(x)
        num = SymbolicUtils.expand(args[1])
        if isadd(num)
            sum(expand_fraction(arg / args[2]) for arg in arguments(num))
        else
            x
        end
    else
        x
    end
end
expand_fraction(x::Num) = Num(expand_fraction(unwrap(x)))

"Apply a function f on every member of a sum or a product"
function _apply_termwise(f, x::BasicSymbolic)
    if isadd(x)
        sum(f(arg) for arg in arguments(x))
    elseif ismul(x)
        prod(f(arg) for arg in arguments(x))
    elseif isdiv(x)
        args = arguments(x)
        _apply_termwise(f, args[1]) / _apply_termwise(f, args[2])
    else
        f(x)
    end
end

simplify_complex(x::Complex) = isequal(x.im, 0) ? x.re : x.re + im * x.im
simplify_complex(x) = x
function simplify_complex(x::BasicSymbolic)
    if isadd(x) || ismul(x) || isdiv(x)
        _apply_termwise(simplify_complex, x)
    else
        # Handle Const-wrapped complex numbers with zero imaginary part
        v = unwrap_const(x)
        if v isa Complex && iszero(imag(v))
            return real(v)
        end
        x
    end
end

"""
$(TYPEDSIGNATURES)

Perform substitutions in `rules` on `x`.
`include_derivatives=true` also substitutes inside derivative arguments using
`Symbolics.substitute_in_deriv`.
"""
Subtype = Union{Num,Equation,BasicSymbolic}

# SymbolicUtils 4 / Symbolics 7: substitute() does not recurse into call arguments
# (e.g. it leaves x(t) unchanged when substituting t => T). Walk the expression and
# rewrite matching nodes ourselves.
function _deep_substitute(e::BasicSymbolic, unwrap_rules::Dict)
    # SU 4's default filter blocks recursion into the arguments of callable-symbolic
    # terms (e.g. `x(t)`). Override with an always-true filter so substitutions like
    # `t => T` reach inside `x(t)`. Compound keys still match before recursion, and
    # SU does not re-walk the replacement, so self-referential rules don't loop.
    return SymbolicUtils.substitute(e, unwrap_rules; filterer=_ -> true)
end
_deep_substitute(x::Num, ur::Dict) = wrap(_deep_substitute(unwrap(x), ur))
function _deep_substitute(eq::Equation, ur::Dict)
    return Equation(
        _deep_substitute(unwrap(eq.lhs), ur), _deep_substitute(unwrap(eq.rhs), ur)
    )
end
_deep_substitute(x, ::Dict) = x

function substitute_all(x::Subtype, rules::Dict; include_derivatives=true)
    unwrap_rules = Dict(unwrap(k) => unwrap(v) for (k, v) in rules)
    result = _deep_substitute(x, unwrap_rules)
    if include_derivatives
        result = Symbolics.substitute_in_deriv(result, rules)
    end
    return result
end

Collections = Union{Dict,Pair,Vector,OrderedDict}
substitute_all(v::AbstractArray, rules) = [substitute_all(x, rules) for x in v]
substitute_all(x::Subtype, rules::Collections) = substitute_all(x, Dict(rules))
function substitute_all(x::Complex{Num}, rules::Collections)
    return substitute_all(x.re, rules) + im * substitute_all(x.im, rules)
end

get_independent(x::Num, t::Num) = wrap(get_independent(unwrap(x), t))
function get_independent(x::Complex{Num}, t::Num)
    return get_independent(x.re, t) + im * get_independent(x.im, t)
end
get_independent(v::Vector{Num}, t::Num) = [get_independent(el, t) for el in v]
get_independent(x, t::Num) = x

function get_independent(x::BasicSymbolic, t::Num)
    if isadd(x)
        sum(get_independent(arg, t) for arg in arguments(x))
    elseif ismul(x)
        prod(get_independent(arg, t) for arg in arguments(x))
    elseif isdiv(x)
        args = arguments(x)
        !is_function(args[2], t) ? get_independent(args[1], t) / args[2] : 0
    elseif ispow(x)
        args = arguments(x)
        !is_function(args[1], t) && !is_function(args[2], t) ? x : 0
    elseif isterm(x) || issym(x)
        !is_function(x, t) ? x : 0
    else
        x
    end
end

"Return all the terms contained in `x`"
get_all_terms(x::Num) = Num.(unique(_get_all_terms(unwrap(Symbolics.expand(x)))))
get_all_terms(x::BasicSymbolic) = unique(_get_all_terms(Symbolics.expand(x)))
function get_all_terms(x::Equation)
    return unique(cat(get_all_terms(Num(x.lhs)), get_all_terms(Num(x.rhs)); dims=1))
end
function _get_all_terms(x::BasicSymbolic)
    if isadd(x)
        vcat([_get_all_terms(arg) for arg in arguments(x)]...)
    elseif ismul(x)
        arguments(x)
    elseif isdiv(x)
        args = arguments(x)
        [_get_all_terms(args[1])..., _get_all_terms(args[2])...]
    else
        [x]
    end
end
_get_all_terms(x) = x

function is_harmonic(x::Num, t::Num)::Bool
    all_terms = get_all_terms(x)
    t_terms = setdiff(all_terms, get_independent(all_terms, t))
    isempty(t_terms) && return true
    trigs = is_trig.(t_terms)

    if !prod(trigs)
        return false
    else
        powers = [max_power(first(arguments(unwrap(term))), t) for term in t_terms[trigs]]
        return all(isone, powers)
    end
end

is_harmonic(x::Equation, t::Num) = is_harmonic(x.lhs, t) && is_harmonic(x.rhs, t)
is_harmonic(x, t) = is_harmonic(Num(x), Num(t))

"Return true if `f` is a function of `var`."
is_function(f, var) = unwrap(var) in get_variables(f)

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
    return unique(_collect_nonpolynomial!(BasicSymbolic[], unwrap(x), unwrapped))
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
_contains_variable(x, vars) = any(v -> isequal(x, v), vars)
function _contains_variable(x::BasicSymbolic, vars)
    any(v -> isequal(x, v), vars) && return true
    SymbolicUtils.iscall(x) || return false
    return any(arg -> _contains_variable(arg, vars), arguments(x))
end

"Walk `x` and push each subexpression which is not polynomial in `vars` onto `terms`."
function _collect_nonpolynomial!(terms, x, vars)
    # a subtree free of `vars` is a coefficient: any operation on it is allowed
    _contains_variable(x, vars) || return terms
    any(v -> isequal(x, v), vars) && return terms # a bare variable

    if isadd(x) || ismul(x)
        for arg in arguments(x)
            _collect_nonpolynomial!(terms, arg, vars)
        end
    elseif ispow(x)
        base, exponent = arguments(x)
        power = unwrap_const(exponent)
        if power isa Real && isinteger(power) && power >= 0
            _collect_nonpolynomial!(terms, base, vars)
        else # negative, fractional or symbolic exponent
            push!(terms, x)
        end
    elseif isdiv(x)
        numerator, denominator = arguments(x)
        _collect_nonpolynomial!(terms, numerator, vars)
        _contains_variable(denominator, vars) && push!(terms, x)
    elseif is_derivative(x)
        _collect_nonpolynomial!(terms, first(arguments(x)), vars)
    else # a call such as sin, exp or log applied to something containing a variable
        push!(terms, x)
    end
    return terms
end

"""
Counts the number of derivatives of a symbolic variable.
"""
function count_derivatives(x::BasicSymbolic)
    (isterm(x) || issym(x)) || error("The input is not a single term or symbol")
    if is_derivative(x)
        # In Symbolics v7, Differential stores the order directly
        op = operation(x)
        return op.order
    else
        return 0
    end
end
count_derivatives(x::Num) = count_derivatives(unwrap(x))
