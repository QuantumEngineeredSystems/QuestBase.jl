"""
    trig_reduce(x)

Simplify trigonometric expressions by converting between exponential and trigonometric forms.
This function performs the following steps:
1. Combines fractions with common denominators
2. Expands all brackets
3. Converts trigonometric functions to exponentials
4. Expands products of exponentials
5. Simplifies exponential products
6. Converts back to trigonometric form
7. Simplifies complex expressions

Returns the simplified expression as a `Num` type.
"""

"""
    is_trig(f::Num)

Check if the given expression `f` is a trigonometric function (sine or cosine).

Returns `true` if `f` is either `sin` or `cos`, `false` otherwise.
"""
function trig_reduce(x::Num)
    combined = expand(add_div(x)) # a/b + c/d = (ad + bc)/bd, then open all brackets
    # Not `expand_all`: its extra `Postwalk(expand_exp_power)` is the most expensive step of
    # the averaging, and `simplify_exp_products` below already normalises `exp(a)^n` at every
    # node it reaches.
    exponential = expand(trig_to_exp(combined)) # expand products of exponentials
    # The exponential form carries complex constants, so the stretch from here to
    # `simplify_complex` is untyped; it is confined to this one expression.
    trigonometric = exp_to_trig(simplify_exp_products(exponential))
    real = real_num(simplify_complex(expand(trigonometric)))
    return reduce_denominator(real) # collapse trig identities left in the denominator
end
trig_reduce(x::SymbolicNode) = trig_reduce(Num(x))

"""
    reduce_denominator(x)

Reduce trigonometric identities left in the denominator of a fraction.

Fraction denominators arise from combining fractions (`add_div` merges `a/b + c/d`
into a single fraction) and from solving linear systems (`rearrange!` divides by the
coefficient determinant). The exponential round-trip in [`trig_reduce`](@ref) only
reaches the numerator, so identities such as the Pythagorean
`cos(ωt)² + sin(ωt)² => 1` survive in the denominator. Reduce the denominator on
its own (it holds no nested fraction, so the recursion terminates) and leave the
numerator untouched. Trig-free denominators are returned unchanged.
"""
reduce_denominator(x::Num) = Num(reduce_denominator(unwrap(x)))
function reduce_denominator(x::SymbolicNode)
    isdiv(x) || return x
    num, den = arguments(x)
    contains_trig(den) || return x
    return num / unwrap(trig_reduce(Num(den)))
end

"Return true if any subterm of `x` is a sin or cos."
contains_trig(x::Num) = contains_trig(unwrap(x))
contains_trig(x) = false
function contains_trig(x::BasicSymbolic)
    is_trig(x) && return true
    SymbolicUtils.iscall(x) || return false
    return any(contains_trig, arguments(x))
end

"Return true if `f` is a sin or cos."
is_trig(f::Num) = is_trig(unwrap(f))
is_trig(f) = false
function is_trig(f::BasicSymbolic)
    f = ispow(f) ? arguments(f)[1] : f
    isterm(f) || return false
    op = operation(f)
    return op === cos || op === sin
end

"""
    add_div(x)

Combine fractions over a common denominator, `a/b + c/d = (ad + bc)/bd`.
"""
add_div(x::Num) = Num(Postwalk(add_with_div)(unwrap(x)))

"""
    trig_to_exp(x::Num)

Convert all trigonometric terms (sin, cos) in expression `x` to their exponential form
using Euler's formula: ``\\exp(ix) = \\cos(x) + i*\\sin(x)``.

Returns the converted expression as a `Num` type.
"""
function trig_to_exp(x::Num)
    all_terms = get_all_terms(x)
    trigs = filter(z -> is_trig(z), all_terms)

    rules = []
    for trig in trigs
        trig_unwrapped = unwrap(trig)
        is_pow = ispow(trig_unwrapped) # trig is either a trig or a power of a trig
        pow_args = is_pow ? arguments(trig_unwrapped) : nothing
        power = is_pow ? unwrap_const(pow_args[2]) : 1  # unwrap Const for numeric power
        base = is_pow ? pow_args[1] : trig_unwrapped
        arg = arguments(base)[1]  # keep as BasicSymbolic for proper symbolic exp(im*arg)
        type = operation(base)

        # Work at BasicSymbolic level to avoid Complex{Num} issues
        if type == cos
            inner = (exp(im * arg) + exp(-im * arg))^power * (1//2)^power
        elseif type == sin
            # sin(x)^n = (im^n) * ((exp(-ix) - exp(ix))/2)^n
            inner = im^power * ((exp(-im * arg) - exp(im * arg)))^power * (1//2)^power
        else
            continue
        end
        term = Num(Symbolics.expand(inner))
        append!(rules, [trig => term])
    end

    result = Symbolics.substitute(x, Dict(rules))
    return convert_to_Num(result)
end
trig_to_exp(x::Complex{Num}) = trig_to_exp(x.re) + im * trig_to_exp(x.im)
convert_to_Num(x::Complex{Num})::Num = Num(first(arguments(unwrap(x.re))))
convert_to_Num(x::Num)::Num = x

"""
    trig_to_exp(x::BasicSymbolic)

Convert all trigonometric terms (sin, cos) in expression `x` to their exponential form
using Euler's formula: ``\\exp(ix) = \\cos(x) + i*\\sin(x)``.
"""
function trig_to_exp(x::BasicSymbolic)
    all_terms = get_all_terms(x)
    trigs = filter(z -> is_trig(z), all_terms)

    rules = []
    for trig in trigs
        is_pow = ispow(trig) # trig is either a trig or a power of a trig
        pow_args = is_pow ? arguments(trig) : nothing
        power = is_pow ? unwrap_const(pow_args[2]) : 1  # unwrap Const for numeric power
        base = is_pow ? pow_args[1] : trig
        arg = arguments(base)[1]
        type = operation(base)

        if type == cos
            term = (exp(im * arg) + exp(-im * arg))^power * (1//2)^power
        elseif type == sin
            term = (1 * im^power) * ((exp(-im * arg) - exp(im * arg)))^power * (1//2)^power
        else
            continue
        end

        append!(rules, [trig => term])
    end
    return Symbolics.substitute(x, Dict(rules))
end

"Convert a single exp(arg) node to trig form with sign normalization.
Returns `nothing` if the node is not an exp call (for use with PassThrough)."
function _exp_to_trig_node(x::BasicSymbolic)
    !(isterm(x) && operation(x) === exp) && return nothing

    arg = first(arguments(x))

    # exp(0) = 1: handle literal zero argument
    arg_val = arg isa BasicSymbolic ? unwrap_const(arg) : arg
    if arg_val isa Number && iszero(arg_val)
        return 1
    end

    trigarg = Symbolics.expand(-im * arg) # the argument of the to-be trig function
    trigarg = simplify_complex(trigarg)

    # After expansion and simplification, check if trigarg is zero (exp(0) = 1)
    ta_val = trigarg isa BasicSymbolic ? unwrap_const(trigarg) : trigarg
    if ta_val isa Number && iszero(ta_val)
        return 1
    end

    # put arguments of trigs into a standard form such that sin(x) = -sin(-x), cos(x) = cos(-x) are recognized
    if isadd(trigarg)
        first_symbol = minimum(
            cat(
                string.(sorted_arguments(trigarg)),
                string.(sorted_arguments(-trigarg));
                dims=1,
            ),
        )

        # put trigarg => -trigarg the lowest alphabetic argument of trigarg is lower than that of -trigarg
        # this is a meaningless key but gives unique signs to all sums
        is_first = minimum(string.(sorted_arguments(trigarg))) == first_symbol
        return if is_first
            cos(-trigarg) - im * sin(-trigarg)
        else
            cos(trigarg) + im * sin(trigarg)
        end
    end
    return if ismul(trigarg) && _has_negative_coefficient(trigarg)
        cos(-trigarg) - im * sin(-trigarg)
    else
        cos(trigarg) + im * sin(trigarg)
    end
end

"""
    exp_to_trig(x::BasicSymbolic)
    exp_to_trig(x)
    exp_to_trig(x::Num)
    exp_to_trig(x::Complex{Num})

Convert exponential expressions to their trigonometric form using
the inverse of Euler's formula:
``\\cos(x) = (\\exp(ix) + \\exp(-ix))/2``
and
``\\sin(x) = (\\exp(ix) - \\exp(-ix))/(2i)``.

Handles various input types including basic symbolic expressions,
complex numbers, and `Num` types. Standardizes the sign of
trigonometric arguments for consistent simplification.
"""
function exp_to_trig(x::BasicSymbolic)
    if isadd(x) || isdiv(x) || ismul(x)
        return _apply_termwise(exp_to_trig, x)
    end
    result = _exp_to_trig_node(x)
    return isnothing(result) ? x : result
end

"Check if a Mul expression has a negative leading coefficient"
function _has_negative_coefficient(x::BasicSymbolic)
    if !ismul(x)
        return false
    end
    # Check the arguments for a negative numeric factor
    # Use unwrap_const to handle Const-wrapped numbers in SymbolicUtils v4
    for arg in arguments(x)
        v = arg isa BasicSymbolic ? unwrap_const(arg) : arg
        if v isa Number && v < 0
            return true
        end
    end
    return false
end

exp_to_trig(x) = x
exp_to_trig(x::Num) = exp_to_trig(unwrap(x))
exp_to_trig(x::Complex{Num}) = exp_to_trig(x.re) + im * exp_to_trig(x.im)
