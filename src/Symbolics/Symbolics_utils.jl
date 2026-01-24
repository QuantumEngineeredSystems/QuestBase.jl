
expand_all(x::Num) = Num(expand_all(x.val))
_apply_termwise(f, x::Num) = wrap(_apply_termwise(f, unwrap(x)))

# Symbolics v7 note: `unwrap` now always returns `BasicSymbolic`; use `Symbolics.value`
# when you want the underlying literal number (via `Const`).
_literal_number(x::Num) = (v = Symbolics.value(x); v isa Number ? v : nothing)
_literal_number(x::BasicSymbolic) =
    (v = SymbolicUtils.unwrap_const(x); v isa Number ? v : nothing)
is_literal_zero(x) = (v = _literal_number(x); v !== nothing && iszero(v))

"Expands using SymbolicUtils.expand and expand_exp_power (changes exp(x)^n to exp(x*n)"
function expand_all(x)
    result = Postwalk(expand_exp_power)(SymbolicUtils.expand(x))
    return isnothing(result) ? x : result
end
expand_all(x::Complex{Num}) = expand_all(x.re) + im * expand_all(x.im)

function expand_fraction(x::BasicSymbolic)
    if isadd(x)
        return _apply_termwise(expand_fraction, x)
    elseif ismul(x)
        return _apply_termwise(expand_fraction, x)
    elseif isdiv(x)
        # Only distribute division over addition in the numerator.
        # In SymbolicUtils v4, `arguments(x.num)` for a multiplication returns factors,
        # so splitting unconditionally would produce incorrect results (e.g. d*(a+b)/c).
        num, den = x.num, x.den
        if isadd(num)
            return sum([expand_fraction(arg) / den for arg in arguments(num)])
        end
        return expand_fraction(num) / expand_fraction(den)
    else
        return x
    end
end
expand_fraction(x::Num) = Num(expand_fraction(x.val))

"Apply a function f on every member of a sum or a product"
function _apply_termwise(f, x::BasicSymbolic)
    if isadd(x)
        return sum([f(arg) for arg in arguments(x)])
    elseif ismul(x)
        return prod([f(arg) for arg in arguments(x)])
    elseif isdiv(x)
        return _apply_termwise(f, x.num) / _apply_termwise(f, x.den)
    else
        return f(x)
    end
end

simplify_complex(x::Complex) = isequal(x.im, 0) ? x.re : x.re + im * x.im
simplify_complex(x) = x
function simplify_complex(x::BasicSymbolic)
    if isadd(x)
        return _apply_termwise(simplify_complex, x)
    elseif ismul(x)
        return _apply_termwise(simplify_complex, x)
    elseif isdiv(x)
        return _apply_termwise(simplify_complex, x)
    else
        return x
    end
end

"""
$(TYPEDSIGNATURES)

Perform substitutions in `rules` on `x`.
`include_derivatives=true` also includes all derivatives of the variables of the keys of `rules`.
"""
Subtype = Union{Num,Equation,BasicSymbolic}
function substitute_all(x::Subtype, rules::Dict; include_derivatives=true)
    if include_derivatives
        drules = Pair[]
        for var in keys(rules)
            if !isa(rules[var], Union{AbstractFloat,Integer})
                pair = Differential(var) => Differential(rules[var])
                push!(drules, pair)
            end
        end
        rules = merge(rules, Dict(drules))
    end
    return substitute(x, rules)
end

Collections = Union{Dict,Pair,Vector,OrderedDict}
substitute_all(v::AbstractArray, rules) = [substitute_all(x, rules) for x in v]
substitute_all(x::Subtype, rules::Collections) = substitute_all(x, Dict(rules))
function substitute_all(x::Complex{Num}, rules::Collections)
    return substitute_all(x.re, rules) + im * substitute_all(x.im, rules)
end

function get_independent(x::Num, t::Num)
    result = get_independent(x.val, t)
    return result isa Complex{Num} ? result : wrap(result)
end
function get_independent(x::Complex{Num}, t::Num)
    return get_independent(x.re, t) + im * get_independent(x.im, t)
end
get_independent(v::Vector{Num}, t::Num) = [get_independent(el, t) for el in v]
get_independent(x, t::Num) = x

function get_independent(x::BasicSymbolic, t::Num)
    if isadd(x)
        return sum([get_independent(arg, t) for arg in arguments(x)])
    elseif ismul(x)
        return prod([get_independent(arg, t) for arg in arguments(x)])
    elseif isdiv(x)
        return !is_function(x.den, t) ? get_independent(x.num, t) / x.den : 0
    elseif ispow(x)
        base, exponent = arguments(x)
        return !is_function(base, t) && !is_function(exponent, t) ? x : 0
    elseif isterm(x)
        return !is_function(x, t) ? x : 0
    elseif issym(x)
        return !is_function(x, t) ? x : 0
    else
        return x
    end
end

"Return all the terms contained in `x`"
get_all_terms(x::Num) = Num.(unique(_get_all_terms(Symbolics.expand(x).val)))
get_all_terms(x::BasicSymbolic) = unique(_get_all_terms(Symbolics.expand(x)))
function get_all_terms(x::Equation)
    return unique(cat(get_all_terms(Num(x.lhs)), get_all_terms(Num(x.rhs)); dims=1))
end
function _get_all_terms(x::BasicSymbolic)
    if isadd(x)
        return vcat([_get_all_terms(term) for term in SymbolicUtils.sorted_arguments(x)]...)
    elseif ismul(x)
        return SymbolicUtils.sorted_arguments(x)
    elseif isdiv(x)
        return [_get_all_terms(x.num)..., _get_all_terms(x.den)...]
    else
        return [x]
    end
end
_get_all_terms(x) = x

function is_harmonic(x::Num, t::Num)::Bool
    all_terms = get_all_terms(x)
    t_terms = setdiff(all_terms, get_independent(all_terms, t))
    isempty(t_terms) && return true
    trigs = is_trig.(t_terms)

        if !all(trigs)
        return false
    else
        powers = [max_power(first(arguments(term.val)), t) for term in t_terms[trigs]]
        return all(isequal(1), powers)
    end
end

is_harmonic(x::Equation, t::Num) = is_harmonic(x.lhs, t) && is_harmonic(x.rhs, t)
is_harmonic(x, t) = is_harmonic(Num(x), Num(t))

"Return true if `f` is a function of `var`."
is_function(f, var) = any(isequal.(get_variables(f), var))

"""
Counts the number of derivatives of a symbolic variable.
"""
function count_derivatives(x::BasicSymbolic)
    if Symbolics.is_derivative(x)
        arg = first(arguments(x))
        (issym(arg) ||
         Symbolics.is_derivative(arg) ||
         (isterm(arg) && issym(operation(arg)))) ||
            error("The input is not a single term or symbol")
        D = operation(x)
        return D.order + count_derivatives(arg)
    end
    (issym(x) || (isterm(x) && issym(operation(x)))) ||
        error("The input is not a single term or symbol")
    return 0
end
count_derivatives(x::Num) = count_derivatives(Symbolics.unwrap(x))
