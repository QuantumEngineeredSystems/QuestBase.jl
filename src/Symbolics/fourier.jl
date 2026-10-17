"""
A variable no expression contains. Projected against it, every trig argument is a phase, so the
time-independent part of the series is the whole expression in product-to-sum form.
"""
const LINEARIZATION_TIME = Symbolics.variable(:questbase_linearization_time)

"""
    trig_reduce(x)

Rewrite `x` with every product and power of `cos` and `sin` expanded into a sum of single
`cos` and `sin` terms, as in `cos(a)^2 => 1/2 + cos(2a)/2`. Identities such as
`cos(a)^2 + sin(a)^2 = 1` therefore cancel, including in denominators, which are reduced on
their own first.

The rewrite is the exact spectral projection behind [`fourier_terms`](@ref): each trig
product is expanded through its exponential spectrum with numeric weights, and the
arguments are written in a canonical sign. Trig arguments may be any linear combination of
symbolic monomials. A trig factor that is not a positive integer power, such as
`1/cos(a)`, is left as it stands.
"""
function trig_reduce(x::Num)
    series = fourier_series(reduce_denominators(x), LINEARIZATION_TIME)
    isempty(series) && return Num(SYMBOLIC_ZERO)
    return only(series).cosine
end
trig_reduce(x::SymbolicNode) = trig_reduce(Num(x))

"`x` with the denominator of every fraction among its summands put through `trig_reduce`."
function reduce_denominators(x::Num)
    summands = summands_of(unwrap(expand(x)))
    any(summand -> contains_trig(last(numerator_denominator(summand))), summands) ||
        return x
    return Num(polynomial(SymbolicNode[reduce_denominator(s) for s in summands]))
end

"""
    reduce_denominator(x)

Reduce trigonometric identities left in the denominator of a fraction.

Fraction denominators arise from solving linear systems (`rearrange!` divides by the
coefficient determinant), and identities such as the Pythagorean `cos(ωt)² + sin(ωt)² => 1`
can survive in them. Reduce the denominator on its own with [`trig_reduce`](@ref) and leave
the numerator untouched. Trig-free denominators are returned unchanged.
"""
reduce_denominator(x::Num) = Num(reduce_denominator(unwrap(x)))
function reduce_denominator(x::SymbolicNode)
    return @match x begin
        BSImpl.Div(; num, den) =>
            contains_trig(den) ? num / unwrap(trig_reduce(Num(den))) : x
        _ => x
    end
end

"Return true if any subterm of `x` is a sin or cos."
contains_trig(x::Num) = contains_trig(unwrap(x))
function contains_trig(x::SymbolicNode)
    is_trig_call(x) && return true
    return @match x begin
        BSImpl.Const(;) || BSImpl.Sym(;) => false
        _ => any(contains_trig, arguments(x))
    end
end

"Return true if `f` is a sin or cos, or a power of one."
is_trig(f::Num) = is_trig(unwrap(f))
is_trig(f::Number) = false
is_trig(f::SymbolicNode) = is_trig_call(first(base_exponent(f)))
