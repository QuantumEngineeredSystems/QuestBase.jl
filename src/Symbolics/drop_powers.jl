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
function drop_powers(expr::BasicSymbolic, vars::AbstractVector, deg::Int)
    deg >= 0 || throw(ArgumentError("deg must be nonnegative"))
    raw_vars = unique(unwrap.(vars))
    all(var -> var isa BasicSymbolic, raw_vars) ||
        throw(ArgumentError("vars must contain symbolic expressions"))
    iszero(deg) && return unwrap(Num(0))
    ϵ = unwrap(Symbolics.variable(gensym(:degree)))
    rules = Dict(var => ϵ * var for var in raw_vars)
    scaled = Symbolics.expand(substitute_all(expr, rules; include_derivatives=false))
    max_deg = max_power(scaled, ϵ)
    removal = Dict(ϵ^d => 0 for d in deg:max_deg)
    truncated = substitute_all(scaled, removal; include_derivatives=false)
    return Symbolics.expand(
        substitute_all(truncated, Dict(ϵ => 1); include_derivatives=false)
    )
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
    terms = get_all_terms(x)
    powers = power_of.(terms, y)
    return maximum(powers)
end

max_power(x::Vector{Num}, y::Num) = maximum(max_power.(x, y))
max_power(x::Complex, y::Num) = maximum(max_power.([x.re, x.im], y))
max_power(x, t) = max_power(wrap(x), wrap(t))

"Return the power of `y` in the term `x`"
function power_of(x::Num, y::Num)
    issym(unwrap(y)) ? nothing : error("power of " * string(y) * " is ambiguous")
    return power_of(unwrap(x), unwrap(y))
end

function power_of(x::BasicSymbolic, y::BasicSymbolic)
    if ispow(x) && issym(y)
        args = arguments(x)
        return isequal(args[1], y) ? unwrap_const(args[2]) : 0
    elseif issym(x) && issym(y)
        return isequal(x, y) ? 1 : 0
    else
        return 0
    end
end

power_of(x, y) = 0
