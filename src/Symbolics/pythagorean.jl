"""
    collapse_pythagorean(x)

Collapse `cos(θ)^2 + sin(θ)^2` to `1` wherever the pair appears in a sum, including with a
shared coefficient, as in `a*cos(θ)^2 + a*sin(θ)^2 => a`.

SymbolicUtils 3 folded this identity away on its own. SymbolicUtils 4 does not, which is
why the determinant-carrying denominators of a symbolic linear solve stopped simplifying
and started growing without bound on a trigonometric ansatz. Applying the identity directly
keeps those expressions small again, at the cost of one bottom-up walk, rather than the full
exponential round trip [`trig_reduce`](@ref) needs.
"""
collapse_pythagorean(x::Num) = wrap(collapse_pythagorean(unwrap(x)))
function collapse_pythagorean(x::Equation)
    return Equation(collapse_pythagorean(x.lhs), collapse_pythagorean(x.rhs))
end
function collapse_pythagorean(x)
    x isa BasicSymbolic || return x
    collapsed = Postwalk(_collapse_pythagorean_node)(x)
    return isnothing(collapsed) ? x : collapsed
end

"`cos`/`sin` and its argument if `x` is one of them, `nothing` otherwise."
function _trig_operation(x)
    (x isa BasicSymbolic && isterm(x)) || return nothing
    op = operation(x)
    (op === cos || op === sin) || return nothing
    return op, first(arguments(x))
end

"Is `x` the literal exponent 2?"
function _is_squared(x)
    value = x isa BasicSymbolic ? unwrap_const(x) : x
    return value isa Number && value == 2
end

"""
Split a summand into `(cos-or-sin, argument, remaining factors)` when it carries a squared
trig factor, `nothing` otherwise.
"""
function _squared_trig_term(summand)
    if summand isa BasicSymbolic && ispow(summand)
        base, exponent = arguments(summand)
        _is_squared(exponent) || return nothing
        trig = _trig_operation(base)
        isnothing(trig) && return nothing
        return trig[1], trig[2], Num(1)
    elseif summand isa BasicSymbolic && ismul(summand)
        trig = nothing
        rest = Num(1)
        for factor in arguments(summand)
            if isnothing(trig) && factor isa BasicSymbolic && ispow(factor)
                base, exponent = arguments(factor)
                if _is_squared(exponent)
                    candidate = _trig_operation(base)
                    if !isnothing(candidate)
                        trig = candidate
                        continue
                    end
                end
            end
            rest *= wrap(factor)
        end
        isnothing(trig) && return nothing
        return trig[1], trig[2], rest
    end
    return nothing
end

"Pair up squared sines and cosines among the summands of a single `Add` node."
function _collapse_pythagorean_node(x)
    (x isa BasicSymbolic && isadd(x)) || return x
    summands = collect(arguments(x))
    length(summands) > 1 || return x
    parts = map(_squared_trig_term, summands)
    any(!isnothing, parts) || return x

    taken = falses(length(summands))
    result = Num(0)
    collapsed_any = false
    for i in eachindex(summands)
        taken[i] && continue
        part = parts[i]
        if !isnothing(part)
            operation_i, argument_i, coefficient_i = part
            partner_operation = operation_i === cos ? sin : cos
            partner = findfirst(eachindex(summands)) do j
                j == i && return false
                taken[j] && return false
                isnothing(parts[j]) && return false
                return parts[j][1] === partner_operation &&
                       isequal(parts[j][2], argument_i) &&
                       isequal(parts[j][3], coefficient_i)
            end
            if !isnothing(partner)
                taken[i] = taken[partner] = true
                result += coefficient_i
                collapsed_any = true
                continue
            end
        end
        taken[i] = true
        result += wrap(summands[i])
    end
    collapsed_any || return x
    return unwrap(result)
end
