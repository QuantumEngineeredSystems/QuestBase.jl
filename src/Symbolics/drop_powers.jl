"""
$(SIGNATURES)
Remove parts of `expr` where the combined power of `vars` is => `deg`.

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
function drop_powers(expr::BasicSymbolic, vars::Vector{<:BasicSymbolic}, deg::Int)
    return unwrap(drop_powers(Num(expr), Num.(vars), deg))
end

# Num fallback: the actual implementation
function drop_powers(expr::Num, vars::Vector{Num}, deg::Int)
    Symbolics.@variables ϵ
    subs_expr = deepcopy(expr)
    rules = Dict([var => ϵ * var for var in unique(vars)])
    subs_expr = Symbolics.expand(
        substitute_all(subs_expr, rules; include_derivatives=false)
    )
    max_deg = max_power(subs_expr, ϵ)
    removal = Dict([ϵ^d => Num(0) for d in deg:max_deg])
    res = substitute_all(
        substitute_all(subs_expr, removal; include_derivatives=false),
        Dict(ϵ => Num(1));
        include_derivatives=false,
    )
    return Symbolics.expand(res)
end

function drop_powers(expr::Vector{Num}, var::Vector{Num}, deg::Int)
    return [drop_powers(x, var, deg) for x in expr]
end

# calls the above for various types of the first argument
function drop_powers(eq::Equation, var::Vector{Num}, deg::Int)
    return drop_powers(eq.lhs, var, deg) .~ drop_powers(eq.rhs, var, deg)
end
function drop_powers(eqs::Vector{Equation}, var::Vector{Num}, deg::Int)
    return [
        Equation(drop_powers(eq.lhs, var, deg), drop_powers(eq.rhs, var, deg)) for eq in eqs
    ]
end
drop_powers(expr::BasicSymbolic, var::BasicSymbolic, deg::Int) = drop_powers(expr, [var], deg)
drop_powers(expr, var::Num, deg::Int) = drop_powers(expr, [var], deg)
drop_powers(x, vars, deg::Int) = drop_powers(Num(x), Num.(vars), deg)

"Return the highest power of `y` occurring in the term `x`."
function max_power(x::BasicSymbolic, y::BasicSymbolic)
    terms = get_all_terms(x)
    powers = power_of.(terms, y)
    return maximum(powers)
end

# Num fallback: unwrap to BasicSymbolic
max_power(x::Num, y::Num) = max_power(unwrap(x), unwrap(y))
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
