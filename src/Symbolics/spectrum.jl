# Fourier coefficients read off exact trigonometric spectra.
#
# After `expand`, every addend is `c * m * Π f_j(a_j)^{p_j}` with `m` free of `t` and each
# `f_j` a `cos` or `sin` whose argument `a_j` is a sum of monomials that all carry `t`.
# Writing `cos(a) = (e^{ia} + e^{-ia})/2` and `sin(a) = (-i e^{ia} + i e^{-ia})/2`, the trig
# part is a finite exponential sum whose weights are plain numbers. Multiplying those small
# tables replaces the symbolic exponential round trip of `trig_reduce`, which expands the
# whole expression twice more and dominates the cost of the averaging.

# A trig argument `Σ c_k m_k` as `m_k => c_k`, with each monomial `m_k` stored as
# `base => exponent`. Two arguments are equal exactly when these dictionaries are.
const _Monomial = Dict{Any,Any}
const _Argument = Dict{_Monomial,Number}
const _Spectrum = Dict{_Argument,Number}

"Split a monomial into its numeric coefficient and `base => exponent` factors."
function _split_monomial(x)
    coefficient = 1
    factors = _Monomial()
    for factor in (ismul(x) ? arguments(x) : (x,))
        value = unwrap_const(factor)
        if value isa Number
            coefficient *= value
            continue
        end
        base, exponent = ispow(factor) ? arguments(factor) : (factor, 1)
        factors[base] = get(factors, base, 0) + unwrap_const(exponent)
    end
    return coefficient, factors
end

"Parse a trig argument linear in monomials that all carry `t`; `nothing` otherwise."
function _trig_argument(a, t)
    argument = _Argument()
    for addend in (isadd(a) ? arguments(a) : (a,))
        coefficient, monomial = _split_monomial(addend)
        # A t-free addend is a phase. Leave it to the general path.
        haskey(monomial, t) || return nothing
        argument[monomial] = get(argument, monomial, 0) + coefficient
    end
    return filter!(!iszero ∘ last, argument)
end

_negate(a::_Argument) = _Argument(k => -v for (k, v) in a)

function _combine(a::_Argument, b::_Argument)
    sum = copy(a)
    for (k, v) in b
        sum[k] = get(sum, k, 0) + v
    end
    return filter!(!iszero ∘ last, sum)
end

function _convolve(s::_Spectrum, f::_Spectrum)
    product = _Spectrum()
    for (a, w) in s, (b, v) in f
        key = _combine(a, b)
        product[key] = get(product, key, 0) + w * v
    end
    return filter!(!iszero ∘ last, product)
end

function _trig_spectrum(f, argument::_Argument)
    w₊, w₋ = f === cos ? (1 // 2, 1 // 2) : (-im // 2, im // 2)
    return _Spectrum(argument => w₊, _negate(argument) => w₋)
end

_is_trig_call(x) = isterm(x) && (operation(x) === cos || operation(x) === sin)

"""
Split an expanded addend into its numeric coefficient, its `t`-free factors and the
spectrum of its trig factors. Return `nothing` if `t` enters other than through
`cos`/`sin` of a phase-free argument raised to a positive integer power.
"""
function _addend_spectrum(x, t)
    coefficient = 1
    rest = Any[]
    spectrum = _Spectrum(_Argument() => 1)
    for factor in (ismul(x) ? arguments(x) : (x,))
        value = unwrap_const(factor)
        if value isa Number
            coefficient *= value
            continue
        end
        base, exponent = ispow(factor) ? arguments(factor) : (factor, 1)
        exponent = unwrap_const(exponent)
        if _is_trig_call(base)
            (exponent isa Integer && exponent > 0) || return nothing
            argument = _trig_argument(unwrap(first(arguments(base))), t)
            isnothing(argument) && return nothing
            factor_spectrum = _trig_spectrum(operation(base), argument)
            for _ in 1:exponent
                spectrum = _convolve(spectrum, factor_spectrum)
            end
        else
            is_function(factor, t) && return nothing
            push!(rest, factor)
        end
    end
    return coefficient, rest, spectrum
end

_real_part(z::Complex) = iszero(imag(z)) ? real(z) : nothing
_real_part(z) = z

"""
    _spectral_terms(x, t, specs)

Coefficients of `f(ω*t)` in `x` for every `(ω, f)` in `specs`, computed from the exact
exponential spectrum of each addend's trig factors. Returns `nothing` when `x` lies
outside that form (a phase in a trig argument, `t` in a denominator or outside a trig
function), so the caller can fall back to the general path.
"""
function _spectral_terms(x, t, specs)
    tu = unwrap(t)
    expression = unwrap(expand(add_div(x)))
    denominator = 1
    if expression isa BasicSymbolic && isdiv(expression)
        expression, denominator = arguments(expression)
        is_function(denominator, t) && return nothing
    end
    targets = map(specs) do (ω, f)
        isequal(ω, 0) && return (_Argument(), nothing)
        argument = _trig_argument(unwrap(expand(ω * t)), tu)
        return isnothing(argument) ? nothing : (argument, f)
    end
    any(isnothing, targets) && return nothing

    terms = [Any[] for _ in specs]
    addends = if expression isa BasicSymbolic && isadd(expression)
        arguments(expression)
    else
        (expression,)
    end
    for addend in addends
        split = _addend_spectrum(addend, tu)
        isnothing(split) && return nothing
        coefficient, rest, spectrum = split
        for (i, (argument, f)) in enumerate(targets)
            weight = if isnothing(f)
                get(spectrum, argument, 0)
            else
                w₊ = get(spectrum, argument, 0)
                w₋ = get(spectrum, _negate(argument), 0)
                f === cos ? w₊ + w₋ : im * (w₊ - w₋)
            end
            iszero(weight) && continue
            # Each addend is a real function of t, so its projections are real.
            weight = _real_part(coefficient * weight)
            isnothing(weight) && return nothing
            term = isempty(rest) ? weight : mul_worker(NUM_VARTYPE, Any[weight; rest])
            # Divide termwise so common factors cancel within each term.
            push!(terms[i], isequal(denominator, 1) ? term : term / denominator)
        end
    end
    return map(terms) do addends
        return expand(wrap(isempty(addends) ? 0 : add_worker(NUM_VARTYPE, addends)))
    end
end
