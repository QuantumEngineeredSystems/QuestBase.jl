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
function trig_reduce(x)
    x = add_div(x) # a/b + c/d = (ad + bc)/bd
    x = expand(x) # open all brackets
    x = trig_to_exp(x)
    # Not `expand_all`: its extra `Postwalk(expand_exp_power)` is the most expensive step of
    # the averaging, and `simplify_exp_products` below already normalises `exp(a)^n` at every
    # node it reaches.
    x = expand(x) # expand products of exponentials
    x = simplify_exp_products(x) # simplify products of exps
    x = exp_to_trig(x)
    x = Num(simplify_complex(expand(x)))
    return reduce_denominator(x) # collapse trig identities left in the denominator
end

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
function reduce_denominator(x)
    u = unwrap(x)
    isdiv(u) || return x
    num, den = arguments(u)
    contains_trig(den) || return x
    return wrap(num) / trig_reduce(wrap(den))
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
    fourier_cos_term(x, ω, t)

Extract the coefficient of cos(ωt) from the expression `x`.
Used in Fourier analysis to find the cosine components of a periodic function.

# Arguments
- `x`: The expression to analyze
- `ω`: The angular frequency
- `t`: The time variable
"""

"""
    add_div(x)

Simplify fractions by combining terms with common denominators.
Transforms expressions of the form a/b + c/d into (ad + bc)/bd.

Returns the simplified fraction as a `Num` type.
"""
function fourier_cos_term(x, ω, t)
    return _fourier_term(x, ω, t, cos)
end

"Simplify fraction a/b + c/d = (ad + bc)/bd"
add_div(x) = wrap(Postwalk(add_with_div)(unwrap(x)))

"""
    fourier_sin_term(x, ω, t)

Extract the coefficient of sin(ωt) from the expression `x`.
Used in Fourier analysis to find the sine components of a periodic function.

# Arguments
- `x`: The expression to analyze
- `ω`: The angular frequency
- `t`: The time variable
"""
function fourier_sin_term(x, ω, t)
    return _fourier_term(x, ω, t, sin)
end

"""
    _fourier_term(x::Equation, ω, t, f)
    _fourier_term(x, ω, t, f)

Internal function to extract Fourier coefficients from expressions.
Handles both equations and regular expressions, returning the coefficient
of the specified trigonometric function f(ωt).

# Arguments
- `x`: The expression or equation to analyze
- `ω`: The angular frequency
- `t`: The time variable
- `f`: The trigonometric function (sin or cos)
"""
function _fourier_term(x::Equation, ω, t, f)
    return Equation(_fourier_term(x.lhs, ω, t, f), _fourier_term(x.rhs, ω, t, f))
end

"Return the coefficient of f(ωt) in `x` where `f` is a cos or sin."
function _fourier_term(x, ω, t, f)
    term = x * f(ω * t)
    term = trig_reduce(term)
    indep = get_independent(term, t)
    # Handle Complex{Num} from get_independent: extract real part
    # (fourier coefficients of real expressions are real)
    if indep isa Complex
        indep = simplify_complex(indep)
        indep = indep isa Complex ? indep.re : indep
    end
    ft = Num(simplify_complex(Symbolics.expand(indep)))
    ft = !isequal(ω, 0) ? 2 * ft : ft # extra factor in case ω = 0 !
    return Symbolics.expand(ft)
end

"""
    fourier_terms(x, t, specs)

Coefficients of `f(ω*t)` in `x` for every `(ω, f)` pair in `specs`, where each `f` is
`cos` or `sin` and `ω == 0` selects the time-independent part.

Equivalent to `[_fourier_term(x, ω, t, f) for (ω, f) in specs]` but expands `x` once
instead of once per pair. Extracting one harmonic with [`_fourier_term`](@ref) costs a
full exponential round trip (`trig_reduce`) of `x * f(ω*t)`, and an ansatz with `n`
harmonics asks for `2n` of them from the same `x`. [`trig_reduce`](@ref) already
linearises `x` into a sum of `coefficient * cos(kωt)` and `coefficient * sin(kωt)`, so
after one call every coefficient can be read straight off that sum.

Falls back to `_fourier_term` per pair whenever the reduced expression is not linear in
the trig atoms, so an input `trig_reduce` cannot fully linearise still gets the right
answer, just at the old cost.
"""
function fourier_terms(x::Equation, t, specs)
    lhs = fourier_terms(x.lhs, t, specs)
    rhs = fourier_terms(x.rhs, t, specs)
    return [Equation(l, r) for (l, r) in zip(lhs, rhs)]
end

function fourier_terms(x, t, specs)
    buckets = _harmonic_buckets(trig_reduce(x), t)
    isnothing(buckets) && return [_fourier_term(x, ω, t, f) for (ω, f) in specs]
    dc, harmonics = buckets
    return [_bucket_term(dc, harmonics, ω, t, f) for (ω, f) in specs]
end

"""
Split a [`trig_reduce`](@ref)d expression into its time-independent part and a
`trig atom => coefficient` mapping. Returns `nothing` if `x` is not linear in trig atoms
whose argument involves `t`, which is the signal to fall back to the per-harmonic path.
"""
function _harmonic_buckets(x, t)
    u = unwrap(x)
    den = Num(1)
    if isdiv(u)
        num, d = arguments(u)
        is_function(d, t) && return nothing # t survives in the denominator
        den = wrap(d)
        u = num
    end
    addends = isadd(u) ? collect(arguments(u)) : [u]
    dc = Num(0)
    harmonics = Dict{Any,Num}()
    for addend in addends
        factors = ismul(addend) ? collect(arguments(addend)) : [addend]
        trigs = filter(z -> is_trig(z) && is_function(z, t), factors)
        if isempty(trigs)
            is_function(addend, t) && return nothing # t-dependent but not through a trig
            dc += wrap(addend) / den
            continue
        end
        length(trigs) == 1 || return nothing # product of harmonics, not linearised
        trig = only(trigs)
        ispow(unwrap(trig)) && return nothing # power of a harmonic, not linearised
        rest = Num(1)
        for factor in factors
            isequal(factor, trig) || (rest *= wrap(factor))
        end
        is_function(rest, t) && return nothing
        harmonics[trig] = get(harmonics, trig, Num(0)) + rest / den
    end
    return dc, harmonics
end

"Read the `f(ω*t)` coefficient out of the buckets built by [`_harmonic_buckets`](@ref)."
function _bucket_term(dc, harmonics, ω, t, f)
    isequal(ω, 0) && return Symbolics.expand(dc)
    for (trig, coefficient) in harmonics
        node = unwrap(trig)
        operation(node) === f || continue
        arg = wrap(first(arguments(node)))
        # trig_reduce normalises signs via cos(-a) = cos(a) and sin(-a) = -sin(a),
        # so the stored atom may carry the negated argument.
        if _vanishes(arg - ω * t)
            return Symbolics.expand(coefficient)
        elseif _vanishes(arg + ω * t)
            return Symbolics.expand(f === cos ? coefficient : -coefficient)
        end
    end
    return Num(0)
end

function _vanishes(x)
    val = unwrap(Symbolics.expand(x))
    val = val isa BasicSymbolic ? unwrap_const(val) : val
    return val isa Number && iszero(val)
end

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
