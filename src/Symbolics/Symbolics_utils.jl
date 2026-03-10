
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
function substitute_all(x::Subtype, rules::Dict; include_derivatives=true)
    result = substitute(x, rules)
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
