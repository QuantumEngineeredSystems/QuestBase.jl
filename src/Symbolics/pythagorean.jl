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

const _TrigSplit = Tuple{Any,Any,Num}

"""
Every way to read a summand as `(cos-or-sin, argument, remaining factors)` by pulling out one
squared trig factor; empty when it carries none.

A summand can hold more than one squared trig factor, as `cos(3θ)^2*sin(θ)^2*a` does, and
only one of its splits pairs up with a partner elsewhere in the sum. Committing to the first
factor found leaves `cos(3θ)^2*sin(θ)^2*a + sin(3θ)^2*sin(θ)^2*a` unpaired, so every split
has to stay a candidate.
"""
function _squared_trig_terms(summand)
    if summand isa BasicSymbolic && ispow(summand)
        base, exponent = arguments(summand)
        _is_squared(exponent) || return _TrigSplit[]
        trig = _trig_operation(base)
        isnothing(trig) && return _TrigSplit[]
        return _TrigSplit[(trig[1], trig[2], Num(1))]
    elseif summand isa BasicSymbolic && ismul(summand)
        factors = collect(arguments(summand))
        splits = _TrigSplit[]
        for (i, factor) in enumerate(factors)
            (factor isa BasicSymbolic && ispow(factor)) || continue
            base, exponent = arguments(factor)
            _is_squared(exponent) || continue
            trig = _trig_operation(base)
            isnothing(trig) && continue

            rest = Num(1)
            for (j, other) in enumerate(factors)
                j == i || (rest *= wrap(other))
            end
            push!(splits, (trig[1], trig[2], rest))
        end
        return splits
    end
    return _TrigSplit[]
end

"Pair up squared sines and cosines among the summands of a single `Add` node."
function _collapse_pythagorean_node(x)
    (x isa BasicSymbolic && isadd(x)) || return x
    summands = collect(arguments(x))
    length(summands) > 1 || return x
    parts = map(_squared_trig_terms, summands)
    any(!isempty, parts) || return x

    taken = falses(length(summands))
    # collected, not accumulated with `+=`: see `_apply_termwise`
    kept = Any[]
    collapsed_any = false
    for i in eachindex(summands)
        taken[i] && continue
        matched = false
        for (operation_i, argument_i, coefficient_i) in parts[i]
            partner_operation = operation_i === cos ? sin : cos
            matches_partner(split) =
                split[1] === partner_operation &&
                isequal(split[2], argument_i) &&
                isequal(split[3], coefficient_i)
            partner = findfirst(eachindex(summands)) do j
                j == i && return false
                taken[j] && return false
                return any(matches_partner, parts[j])
            end
            if !isnothing(partner)
                taken[i] = taken[partner] = true
                push!(kept, unwrap(coefficient_i))
                collapsed_any = true
                matched = true
                break
            end
        end
        matched && continue
        taken[i] = true
        push!(kept, summands[i])
    end
    collapsed_any || return x
    isempty(kept) && return unwrap(Num(0))
    return length(kept) == 1 ? unwrap(Num(only(kept))) : add_worker(NUM_VARTYPE, kept)
end
