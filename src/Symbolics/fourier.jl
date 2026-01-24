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
    x = expand_all(x) # expand products of exponentials
    x = simplify_exp_products(x) # simplify products of exps
    x = exp_to_trig(x)
    x = Num(simplify_complex(expand(x)))
    return x # simplify_fractions(x)# (a*c^2 + b*c)/c^2 = (a*c + b)/c
end

"Return true if `f` is a sin or cos."
is_trig(f::Num) = is_trig(f.val)
is_trig(f) = false
function is_trig(f::BasicSymbolic)
    if ispow(f)
        base, _ = arguments(f)
        f = base
    end
    return isterm(f) && SymbolicUtils.operation(f) ∈ [cos, sin]
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
    ft = Num(simplify_complex(Symbolics.expand(indep)))
    ft = !isequal(ω, 0) ? 2 * ft : ft # extra factor in case ω = 0 !
    return Symbolics.expand(ft)
end

"""
    trig_to_exp(x::Num)

Convert all trigonometric terms (sin, cos) in expression `x` to their exponential form
using Euler's formula: ``\\exp(ix) = \\cos(x) + i*\\sin(x)``.

Returns the converted expression as a `Num` type.
"""
function trig_to_exp(x::Num)
    return Num(trig_to_exp(x.val))
end
trig_to_exp(x::Complex{Num}) = trig_to_exp(x.re) + im * trig_to_exp(x.im)

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
        if is_pow
            base, exponent = arguments(trig)
            power = exponent
            arg = arguments(base)[1]
            type = operation(base)
        else
            power = 1
            arg = arguments(trig)[1]
            type = operation(trig)
        end

        if type == cos
            term = (exp(im * arg) + exp(-im * arg))^power * (1 // 2)^power
        elseif type == sin
            term =
                (1 * im^power) * ((exp(-im * arg) - exp(im * arg)))^power * (1 // 2)^power
        end

        append!(rules, [trig => term])
    end
    return Symbolics.substitute(x, Dict(rules))
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
    elseif isterm(x) && operation(x) === exp
        arg = first(arguments(x))
        trigarg = Symbolics.expand(-im * arg) # the argument of the to-be trig function
        trigarg = simplify_complex(trigarg)

        # put arguments of trigs into a standard form such that sin(x) = -sin(-x), cos(x) = cos(-x) are recognized
        if isadd(trigarg)
            first_symbol = minimum(
                cat(string.(arguments(trigarg)), string.(arguments(-trigarg)); dims=1)
            )

            # put trigarg => -trigarg the lowest alphabetic argument of trigarg is lower than that of -trigarg
            # this is a meaningless key but gives unique signs to all sums
            is_first = minimum(string.(arguments(trigarg))) == first_symbol
            return if is_first
                cos(-trigarg) - im * sin(-trigarg)
            else
                cos(trigarg) + im * sin(trigarg)
            end
        end
        if ismul(trigarg)
            coeff = trigarg.coeff
            if SymbolicUtils.is_literal_number(coeff)
                coeff_val = SymbolicUtils.unwrap_const(coeff)
                if coeff_val isa Real && coeff_val < 0
                    return cos(-trigarg) - im * sin(-trigarg)
                end
            end
        end
        return cos(trigarg) + im * sin(trigarg)
    else
        return x
    end
end

exp_to_trig(x) = x
exp_to_trig(x::Num) = exp_to_trig(x.val)
exp_to_trig(x::Complex{Num}) = exp_to_trig(x.re) + im * exp_to_trig(x.im)
