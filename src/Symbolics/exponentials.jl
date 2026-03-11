"Returns true if expr is an exponential"
isexp(expr::BasicSymbolic) = @match expr begin
    BSImpl.Term(; f) => f === exp
    _ => false
end
isexp(expr) = false

"Expand powers of exponential such that exp(x)^n => exp(x*n)"
function expand_exp_power(expr::BasicSymbolic)
    @match expr begin
        BSImpl.AddMul(; variant) => if variant == AddMulVariant.ADD
            sum(expand_exp_power(arg) for arg in arguments(expr))
        else  # MUL
            prod(expand_exp_power(arg) for arg in arguments(expr))
        end
        BSImpl.Term(; f, args) => if f === (^) && isexp(args[1])
            exp(arguments(args[1])[1] * args[2])
        else
            expr
        end
        _ => expr
    end
end
expand_exp_power(expr::Num) = expand_exp_power(unwrap(expr))
expand_exp_power(expr) = expr

"Simplify products of exponentials such that exp(a)*exp(b) => exp(a+b).
Also expands exp(x)^n => exp(x*n) and simplifies exp(0) => 1."
function simplify_exp_products(expr::BasicSymbolic)
    @match expr begin
        BSImpl.AddMul(; variant) => if variant == AddMulVariant.ADD
            _apply_termwise(simplify_exp_products, expr)
        else  # MUL
            _simplify_exp_products_mul(expr)
        end
        BSImpl.Div() => _apply_termwise(simplify_exp_products, expr)
        BSImpl.Term(; f, args) => if f === (^) && isexp(args[1])
            # exp(x)^n => exp(x*n)
            exp(arguments(args[1])[1] * args[2])
        else
            expr
        end
        _ => expr
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
