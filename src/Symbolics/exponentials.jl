"Returns true if expr is an exponential"
isexp(expr) = isterm(expr) && operation(expr) === exp

"""
The argument `a` of an exponential factor, reading `exp(a)` and `exp(a)^n` alike (the latter
as `a*n`); `nothing` otherwise.

A `Mul` stores repeated factors as a power, so `exp(a)*exp(a)` is held as `exp(a)^2` and a
plain [`isexp`](@ref) test misses it.
"""
function _exp_argument(expr)
    isexp(expr) && return arguments(expr)[1]
    if expr isa BasicSymbolic && ispow(expr)
        base, exponent = arguments(expr)
        isexp(base) && return arguments(base)[1] * exponent
    end
    return nothing
end

"Expand powers of exponential such that exp(x)^n => exp(x*n)"
function expand_exp_power(expr::BasicSymbolic)
    if isadd(expr)
        add_worker(vartype(typeof(expr)), map(expand_exp_power, arguments(expr)))
    elseif ismul(expr)
        mul_worker(vartype(typeof(expr)), map(expand_exp_power, arguments(expr)))
    elseif ispow(expr) && isexp(arguments(expr)[1])
        exp(arguments(arguments(expr)[1])[1] * arguments(expr)[2])
    else
        expr
    end
end
expand_exp_power(expr::Num) = expand_exp_power(unwrap(expr))
expand_exp_power(expr) = expr

"Simplify products of exponentials such that exp(a)*exp(b) => exp(a+b).
Also expands exp(x)^n => exp(x*n) and simplifies exp(0) => 1."
function simplify_exp_products(expr::BasicSymbolic)
    if isadd(expr)
        _apply_termwise(simplify_exp_products, expr)
    elseif isdiv(expr)
        _apply_termwise(simplify_exp_products, expr)
    elseif ismul(expr)
        _simplify_exp_products_mul(expr)
    elseif ispow(expr) && isexp(arguments(expr)[1])
        # exp(x)^n => exp(x*n), which is otherwise left unexpanded here
        exp(arguments(arguments(expr)[1])[1] * arguments(expr)[2])
    else
        expr
    end
end

function _simplify_exp_products_mul(expr)
    args = arguments(expr)
    exp_args = map(_exp_argument, args)
    ind = findall(!isnothing, exp_args)
    rest_ind = setdiff(1:length(args), ind)
    vtype = vartype(typeof(expr))
    rest = isempty(rest_ind) ? 1 : mul_worker(vtype, args[rest_ind])
    total = isempty(ind) ? 0 : add_worker(vtype, exp_args[ind])
    if is_literal_number(total)
        iszero(unwrap_const(total)) && return rest
    end
    return isempty(ind) ? rest : rest * exp(total)
end

simplify_exp_products(x::Num) = simplify_exp_products(unwrap(x))
function simplify_exp_products(x::Complex{Num})
    return Complex{Num}(
        simplify_exp_products(unwrap(x.re)), simplify_exp_products(unwrap(x.im))
    )
end
simplify_exp_products(x) = x
