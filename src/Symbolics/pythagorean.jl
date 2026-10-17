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
collapse_pythagorean(x::Num) = Num(collapse_pythagorean(unwrap(x)))
function collapse_pythagorean(x::Equation)
    return Equation(collapse_pythagorean(x.lhs), collapse_pythagorean(x.rhs))
end
collapse_pythagorean(x::SymbolicNode) = Postwalk(collapse_pythagorean_sum)(x)

"""
One way to read summand `summand` of a sum as `coefficient * f(argument)^2`, with `f` the
cosine when `iscos` and the sine otherwise.
"""
struct SquaredTrig
    summand::Int
    iscos::Bool
    argument::SymbolicNode
    coefficient::SymbolicNode
end

"Key under which a [`SquaredTrig`](@ref) waits for its partner."
squared_trig_key(split::SquaredTrig) = (split.iscos, split.argument, split.coefficient)
"Key of the partner that completes `split` to `coefficient * (cos^2 + sin^2)`."
partner_key(split::SquaredTrig) = (!split.iscos, split.argument, split.coefficient)

"Is `x` a cosine or sine raised to the literal power 2?"
function is_squared_trig(x::SymbolicNode)
    ispow(x) || return false
    base, exponent = arguments(x)
    # compared as a node: `unwrap_const` is untyped, and `isequal` matches 2, 2.0 and 2//1 alike
    isequal(exponent, SYMBOLIC_TWO) || return false
    isterm(base) || return false
    op = operation(base)
    return op === cos || op === sin
end

"Append the reading of a bare `f(θ)^2` summand, with coefficient one."
function push_squared_trig!(splits::Vector{SquaredTrig}, index::Int, power::SymbolicNode)
    base = first(arguments(power))
    push!(
        splits,
        SquaredTrig(index, operation(base) === cos, first(arguments(base)), SYMBOLIC_ONE),
    )
    return splits
end

"""
Append every way to read summand `index` as `coefficient * f(θ)^2` by pulling out one squared
trig factor; nothing is appended when it carries none.

A summand can hold more than one squared trig factor, as `cos(3θ)^2*sin(θ)^2*a` does, and
only one of its splits pairs up with a partner elsewhere in the sum. Committing to the first
factor found leaves `cos(3θ)^2*sin(θ)^2*a + sin(3θ)^2*sin(θ)^2*a` unpaired, so every split
has to stay a candidate.
"""
function append_squared_trigs!(
    splits::Vector{SquaredTrig}, index::Int, summand::SymbolicNode
)
    is_squared_trig(summand) && return push_squared_trig!(splits, index, summand)
    ismul(summand) || return splits
    factors = arguments(summand)
    for (i, factor) in enumerate(factors)
        is_squared_trig(factor) || continue
        rest = SymbolicNode[factors[j] for j in eachindex(factors) if j != i]
        base = first(arguments(factor))
        coefficient = length(rest) == 1 ? only(rest) : mul_worker(NUM_VARTYPE, rest)
        push!(
            splits,
            SquaredTrig(
                index, operation(base) === cos, first(arguments(base)), coefficient
            ),
        )
    end
    return splits
end

"""
Pair up squared sines and cosines among the summands of a single sum node.

Each split is filed under `(iscos, argument, coefficient)`, so finding the partner of a split
is one hash lookup rather than a scan over the other summands. SymbolicUtils hash-conses its
expressions, so hashing the argument and coefficient is cheap. Summands are visited in order
and each pairs with the earliest summand still free, which keeps the result deterministic.
"""
function collapse_pythagorean_sum(x::SymbolicNode)
    isadd(x) || return x
    summands = arguments(x)
    splits = SquaredTrig[]
    for (index, summand) in enumerate(summands)
        append_squared_trigs!(splits, index, summand)
    end
    length(splits) > 1 || return x

    candidates = Dict{Tuple{Bool,SymbolicNode,SymbolicNode},Vector{Int}}()
    for split in splits
        push!(get!(Vector{Int}, candidates, squared_trig_key(split)), split.summand)
    end

    taken = falses(length(summands))
    # collected, not accumulated with `+=`: see `_apply_termwise`
    kept = SymbolicNode[]
    collapsed_any = false
    next_split = 1
    for index in eachindex(summands)
        # splits are appended in summand order, so those of `index` are contiguous
        first_split = next_split
        while next_split <= length(splits) && splits[next_split].summand == index
            next_split += 1
        end
        taken[index] && continue
        taken[index] = true
        partner = 0
        coefficient = SYMBOLIC_ONE
        for split in view(splits, first_split:(next_split - 1))
            haskey(candidates, partner_key(split)) || continue
            for j in candidates[partner_key(split)]
                taken[j] && continue
                partner = j
                break
            end
            partner == 0 && continue
            coefficient = split.coefficient
            break
        end
        if partner == 0
            push!(kept, summands[index])
        else
            taken[partner] = true
            push!(kept, coefficient)
            collapsed_any = true
        end
    end
    collapsed_any || return x
    return length(kept) == 1 ? only(kept) : add_worker(NUM_VARTYPE, kept)
end
