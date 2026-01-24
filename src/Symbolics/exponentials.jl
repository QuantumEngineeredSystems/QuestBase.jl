expand_exp_power(expr::Num) = expand_exp_power(expr.val)
simplify_exp_products(x::Num) = simplify_exp_products(x.val)

"Returns true if expr is an exponential"
isexp(expr) = isterm(expr) && operation(expr) === exp

"Expand powers of exponential such that exp(x)^n => exp(x*n) "
function expand_exp_power(expr::BasicSymbolic)
    if isadd(expr)
        return sum([expand_exp_power(arg) for arg in arguments(expr)])
    elseif ismul(expr)
        return prod([expand_exp_power(arg) for arg in arguments(expr)])
    else
        if ispow(expr)
            base, exponent = arguments(expr)
            if isexp(base)
                return exp(arguments(base)[1] * exponent)
            end
        end
        return expr
    end
end
expand_exp_power(expr) = expr

"Simplify products of exponentials such that exp(a)*exp(b) => exp(a+b)
This is included in SymbolicUtils as of 17.0 but the method here avoid other simplify calls"
function simplify_exp_products(expr::BasicSymbolic)
    if isadd(expr)
        return _apply_termwise(simplify_exp_products, expr)
    elseif isdiv(expr)
        return _apply_termwise(simplify_exp_products, expr)
    elseif ismul(expr)
        return simplify_exp_products_mul(expr)
    else
        return expr
    end
end
function simplify_exp_products(x::Complex{Num})
    return Complex{Num}(simplify_exp_products(x.re.val), simplify_exp_products(x.im.val))
end
function simplify_exp_products_mul(expr)
    ind = findall(x -> isexp(x), arguments(expr))
    rest_ind = setdiff(1:length(arguments(expr)), ind)
    rest = isempty(rest_ind) ? 1 : prod(arguments(expr)[rest_ind])
    total = isempty(ind) ? 0 : sum(getindex.(arguments.(arguments(expr)[ind]), 1))
    if SymbolicUtils.is_literal_number(total)
        return iszero(SymbolicUtils.unwrap_const(total)) ? rest : rest * exp(total)
    end
    return rest * exp(total)
end
simplify_exp_products(x) = x
