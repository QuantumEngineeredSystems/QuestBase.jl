expand_exp_power(expr::Num) = expand_exp_power(unwrap(expr))
simplify_exp_products(x::Num) = simplify_exp_products(unwrap(x))

"Returns true if expr is an exponential"
isexp(expr) = isterm(expr) && operation(expr) === exp

"Expand powers of exponential such that exp(x)^n => exp(x*n) "
function expand_exp_power(expr::BasicSymbolic)
    if isadd(expr)
        sum(expand_exp_power(arg) for arg in arguments(expr))
    elseif ismul(expr)
        prod(expand_exp_power(arg) for arg in arguments(expr))
    else
        if ispow(expr) && isexp(arguments(expr)[1])
            exp(arguments(arguments(expr)[1])[1] * arguments(expr)[2])
        else
            expr
        end
    end
end
expand_exp_power(expr) = expr

"Simplify products of exponentials such that exp(a)*exp(b) => exp(a+b)
This is included in SymbolicUtils as of 17.0 but the method here avoid other simplify calls"
function simplify_exp_products(expr::BasicSymbolic)
    if isadd(expr)
        _apply_termwise(simplify_exp_products, expr)
    elseif isdiv(expr)
        _apply_termwise(simplify_exp_products, expr)
    elseif ismul(expr)
        simplify_exp_products_mul(expr)
    else
        expr
    end
end
function simplify_exp_products(x::Complex{Num})
    return Complex{Num}(
        simplify_exp_products(unwrap(x.re)), simplify_exp_products(unwrap(x.im))
    )
end
function simplify_exp_products_mul(expr)
    ind = findall(x -> isexp(x), arguments(expr))
    rest_ind = setdiff(1:length(arguments(expr)), ind)
    rest = isempty(rest_ind) ? 1 : prod(arguments(expr)[rest_ind])
    total = isempty(ind) ? 0 : sum(getindex.(arguments.(arguments(expr)[ind]), 1))
    if is_literal_number(total)
        (iszero(unwrap_const(total)) && return rest)
    else
        return rest * exp(total)
    end
end
simplify_exp_products(x) = x
