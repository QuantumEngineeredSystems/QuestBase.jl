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
    x = _trig_expand_products(x)
    x = Num(simplify_complex(expand(x)))
    return x # simplify_fractions(x)# (a*c^2 + b*c)/c^2 = (a*c + b)/c
end

function _is_sin_cos(ex::BasicSymbolic)
    return isterm(ex) && (operation(ex) === sin || operation(ex) === cos)
end

function _trig_mul_to_sum(a::BasicSymbolic, b::BasicSymbolic)
    op1, op2 = operation(a), operation(b)
    x = first(arguments(a))
    y = first(arguments(b))
    if op1 === cos && op2 === cos
        return (cos(x - y) + cos(x + y)) / 2
    elseif op1 === sin && op2 === sin
        return (cos(x - y) - cos(x + y)) / 2
    elseif op1 === sin && op2 === cos
        return (sin(x + y) + sin(x - y)) / 2
    elseif op1 === cos && op2 === sin
        return (sin(x + y) - sin(x - y)) / 2
    end
    return nothing
end

function _trig_expand_products(x::BasicSymbolic)
    # Expand trig products/powers into sums so `get_independent` can isolate constants.
    y = Postwalk(ex -> begin
        if ispow(ex)
            base, exponent = arguments(ex)
            exp_val = SymbolicUtils.unwrap_const(exponent)
            if exp_val isa Integer && exp_val == 2 && _is_sin_cos(base)
                arg = first(arguments(base))
                if operation(base) === cos
                    return (1 + cos(2 * arg)) / 2
                else
                    return (1 - cos(2 * arg)) / 2
                end
            end
        elseif ismul(ex)
            # In SymbolicUtils v4, `arguments(ismul(...))` includes the numeric coefficient
            # even though `ex.coeff` also stores it. Avoid double-counting it.
            factors = BasicSymbolic[
                f for f in arguments(ex) if !(SymbolicUtils.unwrap_const(f) isa Number)
            ]
            trig_idx = findall(_is_sin_cos, factors)
            if length(trig_idx) >= 2
                i, j = trig_idx[1], trig_idx[2]
                repl = _trig_mul_to_sum(factors[i], factors[j])
                if repl !== nothing
                    others = BasicSymbolic[]
                    for (k, f) in pairs(factors)
                        (k == i || k == j) && continue
                        push!(others, f)
                    end
                    coeff = ex.coeff
                    return coeff * prod(others; init=1) * repl
                end
            end
        end
        return ex
    end)(x)
    return SymbolicUtils.expand(y)
end
_trig_expand_products(x::Num) = wrap(_trig_expand_products(unwrap(x)))
_trig_expand_products(x) = x

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

_real_if_complex(v) = v isa Complex{Num} ? v.re : v

_canonicalize(ft) = _real_if_complex(Symbolics.simplify(Symbolics.expand(ft)))

function _cleanup_fourier_term(ft)
    ft = _real_if_complex(ft)
    ft = _strip_real_imag(ft)
    ft = _canonicalize(ft)
    ft = _simplify_trig_zero(Num(ft))
    ft = simplify_complex(Symbolics.expand(ft))
    ft = _real_if_complex(ft)
    ft = _strip_real_imag(ft)
    ft = _normalize_trig_signs(unwrap(ft))
    ft = _strip_zero_imag_literals(wrap(ft))
    ft = _canonicalize(ft)
    ft = _real_if_complex(ft)
    return ft
end

"Return the coefficient of f(ωt) in `x` where `f` is a cos or sin."
function _fourier_term(x, ω, t, f)
    term = x * f(ω * t)
    term = trig_reduce(term)
    indep = get_independent(term, t)
    ft = simplify_complex(Symbolics.expand(indep))
    ft = !isequal(ω, 0) ? 2 * ft : ft # extra factor in case ω = 0 !
    ft = _cleanup_fourier_term(ft)
    return ft isa Num ? ft : Num(ft)
end

function _strip_real_imag(x::Complex{Num})
    return _strip_real_imag(x.re) + im * _strip_real_imag(x.im)
end

function _strip_zero_imag_literals(x::BasicSymbolic)
    return Postwalk(ex -> begin
        v = SymbolicUtils.unwrap_const(ex)
        if v isa Complex && iszero(imag(v))
            return real(v)
        end
        return ex
    end)(x)
end

function _strip_zero_imag_literals(x::Num)
    return wrap(_strip_zero_imag_literals(unwrap(x)))
end

function _strip_zero_imag_literals(x::Complex{Num})
    return _strip_zero_imag_literals(x.re) + im * _strip_zero_imag_literals(x.im)
end

function _strip_real_imag(x::BasicSymbolic)
    function _real_of(ex::BasicSymbolic)
        if isadd(ex)
            return sum(_real_of.(arguments(ex)))
        elseif ismul(ex)
            coeff_val = SymbolicUtils.unwrap_const(ex.coeff)
            if coeff_val isa Number
                r = real(coeff_val)
                rest = prod(
                    (f for f in arguments(ex) if !(SymbolicUtils.unwrap_const(f) isa Number));
                    init=1,
                )
                return iszero(r) ? 0 : r * rest
            end
            return ex
        elseif isdiv(ex)
            return _real_of(ex.num) / ex.den
        else
            v = SymbolicUtils.unwrap_const(ex)
            return v isa Number ? real(v) : ex
        end
    end

    function _imag_of(ex::BasicSymbolic)
        if isadd(ex)
            return sum(_imag_of.(arguments(ex)))
        elseif ismul(ex)
            coeff_val = SymbolicUtils.unwrap_const(ex.coeff)
            if coeff_val isa Number
                i = imag(coeff_val)
                rest = prod(
                    (f for f in arguments(ex) if !(SymbolicUtils.unwrap_const(f) isa Number));
                    init=1,
                )
                return iszero(i) ? 0 : i * rest
            end
            return 0
        elseif isdiv(ex)
            return _imag_of(ex.num) / ex.den
        else
            v = SymbolicUtils.unwrap_const(ex)
            return v isa Number ? imag(v) : 0
        end
    end

    return Postwalk(ex -> begin
        if isterm(ex)
            op = operation(ex)
            if op === real || nameof(op) == :real
                return _real_of(first(arguments(ex)))
            elseif op === imag || nameof(op) == :imag
                return _imag_of(first(arguments(ex)))
            end
        end
        return ex
    end)(x)
end

function _strip_real_imag(x::Num)
    return wrap(_strip_real_imag(unwrap(x)))
end

function _simplify_trig_zero(x::BasicSymbolic)
    return Postwalk(ex -> begin
        if isterm(ex)
            op = operation(ex)
            if op === sin || op === cos
                arg = first(arguments(ex))
                arg_val = SymbolicUtils.unwrap_const(arg)
                if arg_val isa Number && iszero(arg_val)
                    return op === sin ? 0 : 1
                end
            end
        end
        return ex
    end)(x)
end

function _simplify_trig_zero(x::Num)
    return wrap(_simplify_trig_zero(unwrap(x)))
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
        return Symbolics.simplify(_normalize_trig_signs(_apply_termwise(exp_to_trig, x)))
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
            coeff_val = SymbolicUtils.unwrap_const(coeff)
            if coeff_val isa Real && coeff_val < 0
                return cos(-trigarg) - im * sin(-trigarg)
            end
        end
        return _normalize_trig_signs(cos(trigarg) + im * sin(trigarg))
    else
        return x
    end
end

function _normalize_trig_signs(x::BasicSymbolic)
    if isadd(x) || ismul(x)
        args = SymbolicUtils.sorted_arguments(x)
        return (isadd(x) ? sum : prod)(_normalize_trig_signs.(args))
    elseif isdiv(x)
        return _normalize_trig_signs(x.num) / _normalize_trig_signs(x.den)
    elseif isterm(x)
        op = operation(x)
        if op === real || op === imag || nameof(op) == :real || nameof(op) == :imag
            arg = first(arguments(x))
            return op(_normalize_trig_signs(arg))
        elseif op === sin || op === cos
            arg = first(arguments(x))
            if SymbolicUtils.isnegative(arg)
                new_arg = -arg
                return op === sin ? -sin(new_arg) : cos(new_arg)
            elseif ismul(arg)
                coeff_val = SymbolicUtils.unwrap_const(arg.coeff)
                # Handle coefficients that are complex with zero imaginary part, e.g. (-2 + 0im)θ.
                if coeff_val isa Number && isreal(coeff_val) && real(coeff_val) < 0
                    new_arg = -arg
                    return op === sin ? -sin(new_arg) : cos(new_arg)
                end
            end
        end
    end
    return x
end

exp_to_trig(x) = x
exp_to_trig(x::Num) = exp_to_trig(x.val)
exp_to_trig(x::Complex{Num}) = exp_to_trig(x.re) + im * exp_to_trig(x.im)
