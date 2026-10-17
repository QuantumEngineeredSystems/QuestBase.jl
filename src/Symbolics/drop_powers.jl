"""
$(SIGNATURES)
Remove terms of `expr` whose combined degree in `vars` is at least `deg`.

Accept raw `BasicSymbolic` expressions or `Num` wrappers, with a symbolic variable or
a vector of variables. Raw symbolic inputs stay raw; `Num` inputs return `Num`.
Vectors and equations are truncated componentwise.
`deg` must be nonnegative; degree zero removes every term.

# Example
```julia-repl
julia> @variables x,y;
julia>drop_powers((x+y)^2, x, 2)
y^2 + 2*x*y
julia>drop_powers((x+y)^2, [x,y], 2)
0
julia>drop_powers((x+y)^2 + (x+y)^3, [x,y], 3)
x^2 + y^2 + 2*x*y
```
"""
function drop_powers(expr::SymbolicNode, vars::AbstractVector, deg::Int)
    deg >= 0 || throw(ArgumentError("deg must be nonnegative"))
    all(var -> unwrap(var) isa BasicSymbolic, vars) ||
        throw(ArgumentError("vars must contain symbolic expressions"))
    iszero(deg) && return SYMBOLIC_ZERO
    raw_vars = unique(SymbolicNode[unwrap(var) for var in vars])
    terms = expanded_terms(expr)
    return expanded_expression(filter!(term -> !is_dropped(term, raw_vars, deg), terms))
end

"""
Whether the expanded `term` has combined degree at least `deg` in `vars`. A fraction's
degree is that of its numerator less that of its denominator. A term holding a power of a
variable that is not a number has no degree and is kept.
"""
function is_dropped(term::ExpandedTerm, vars::Vector{SymbolicNode}, deg::Int)
    num_degree, num_determined = monomial_degree(term.factors, vars)
    den_degree, den_determined = monomial_degree(factors_of(term.denominator), vars)
    degree = num_degree - den_degree
    return num_determined && den_determined && isinteger(degree) && degree >= deg
end

"""
The combined degree in `vars` of `factors`, and whether every power of a variable among
them is a number. A factor that is not a power of a variable, such as `sin(x)` or a sum,
has degree zero.
"""
function monomial_degree(factors, vars::Vector{SymbolicNode})
    degree = 0//1
    for factor in factors
        base, exponent = base_exponent(factor)
        any(var -> isequal(base, var), vars) || continue
        isconst(exponent) || return degree, false
        degree += with_const(Rational{Int}, exponent)
    end
    return degree, true
end

function drop_powers(expr::Num, vars::AbstractVector, deg::Int)
    return Num(drop_powers(unwrap(expr), vars, deg))
end
drop_powers(expr::Real, vars::AbstractVector, deg::Int) = drop_powers(Num(expr), vars, deg)
function drop_powers(expr::AbstractVector, vars::AbstractVector, deg::Int)
    return [drop_powers(term, vars, deg) for term in expr]
end
function drop_powers(eq::Equation, vars::AbstractVector, deg::Int)
    return Equation(drop_powers(eq.lhs, vars, deg), drop_powers(eq.rhs, vars, deg))
end
drop_powers(expr, var::Union{Num,BasicSymbolic}, deg::Int) = drop_powers(expr, [var], deg)

"Return the highest power of `y` occurring in the term `x`."
function max_power(x::Num, y::Num)
    return maximum(term -> power_of(term, y), get_all_terms(x))
end

max_power(x::Vector{Num}, y::Num) = maximum(term -> max_power(term, y), x)
max_power(x::Complex{Num}, y::Num) = max(max_power(x.re, y), max_power(x.im, y))
max_power(x::SymbolicNode, y::SymbolicNode) = max_power(Num(x), Num(y))

"Return the power of the symbol `y` in the term `x`."
function power_of(x::Num, y::Num)
    issym(unwrap(y)) || throw(ArgumentError("The power of $y is ambiguous."))
    return power_of(unwrap(x), unwrap(y))
end

function power_of(x::SymbolicNode, y::SymbolicNode)
    issym(y) || return 0//1
    base, exponent = base_exponent(x)
    isequal(base, y) || return 0//1
    isconst(exponent) || throw(ArgumentError("The power of $y in $x is not a number."))
    return with_const(Rational{Int}, exponent)
end
