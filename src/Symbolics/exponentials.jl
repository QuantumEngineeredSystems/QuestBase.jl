"Returns true if expr is an exponential"
isexp(expr) = isterm(expr) && operation(expr) === exp

"Expand powers of exponential such that exp(x)^n => exp(x*n)"
function expand_exp_power(expr::BasicSymbolic)
    if isadd(expr)
        sum(expand_exp_power(arg) for arg in arguments(expr))
    elseif ismul(expr)
        prod(expand_exp_power(arg) for arg in arguments(expr))
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
        # exp(x)^n => exp(x*n) — fixes bug where exp powers were left unexpanded
        exp(arguments(arguments(expr)[1])[1] * arguments(expr)[2])
    else
        expr
    end
end

function _simplify_exp_products_mul(expr)
    args = arguments(expr)
    ind = findall(isexp, args)
    rest_ind = setdiff(1:length(args), ind)
    rest = isempty(rest_ind) ? 1 : prod(args[rest_ind])
    total = isempty(ind) ? 0 : sum(arguments(args[i])[1] for i in ind)
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
