# Distributive expansion into plain data.
#
# `SymbolicUtils.expand` builds a hash-consed node for every partial product, and those
# constructions, not the arithmetic, are what it spends its time on. Here products and
# powers of sums are distributed into `ExpandedTerm`s, which hold their factors in a vector,
# so no node is built until a caller asks for one. The Fourier projection and `drop_powers`
# read the terms directly and build nodes only for the terms they return.

"""
One term `Π factors / denominator` of an expanded expression. No factor is a sum or a
fraction; a numeric constant may appear among the factors more than once.
"""
struct ExpandedTerm
    factors::Vector{SymbolicNode}
    denominator::SymbolicNode
end

ExpandedTerm(factor::SymbolicNode) = ExpandedTerm(SymbolicNode[factor], SYMBOLIC_ONE)

"The product of two denominators, either of which may be one."
function denominator_product(a::SymbolicNode, b::SymbolicNode)
    isequal(a, SYMBOLIC_ONE) && return b
    isequal(b, SYMBOLIC_ONE) && return a
    return mul_worker(NUM_VARTYPE, SymbolicNode[a, b])
end

function Base.:*(a::ExpandedTerm, b::ExpandedTerm)
    return ExpandedTerm(
        vcat(a.factors, b.factors), denominator_product(a.denominator, b.denominator)
    )
end

"Every product of a term of `a` with a term of `b`."
function distribute(a::Vector{ExpandedTerm}, b::Vector{ExpandedTerm})
    products = Vector{ExpandedTerm}(undef, length(a) * length(b))
    k = 0
    for x in a, y in b
        products[k += 1] = x * y
    end
    return products
end

"Return true if the constant `value` is an integer power that expansion distributes."
is_distributed_power(value::Integer) = value > 1
is_distributed_power(value::Number) = false

"Return true if the constant `value` is a negative integer, a power that is a denominator."
is_inverse_power(value::Integer) = value < 0
is_inverse_power(value::Number) = false

"The integer power `value` as a count of repeated factors."
repetitions(value::Integer) = Int(abs(value))
repetitions(value::Number) = throw(ArgumentError("$value is not an integer power."))

"The number `value` as a symbolic constant."
constant_node(value::Number) = Const{NUM_VARTYPE}(value)

"""
    expanded_terms(x)

The terms of `x` with every product and positive integer power of a sum distributed, and
every fraction split over the terms of its numerator. Denominators are expanded on their
own and kept whole. Arguments of function calls are left as they stand.

Sums and products are read from their coefficient dictionaries rather than through
`arguments`, which would build a node for every summand and power. A positive integer power
of an atom enters as that many copies of the atom.
"""
function expanded_terms(x::SymbolicNode)
    return @match x begin
        BSImpl.AddMul(; coeff, dict, variant) => if variant === AddMulVariant.ADD
            terms = ExpandedTerm[]
            with_number(iszero, coeff) ||
                push!(terms, ExpandedTerm(with_number(constant_node, coeff)))
            for (summand, weight) in dict
                append_weighted!(terms, expanded_terms(summand), weight)
            end
            terms
        else
            terms = [ExpandedTerm(SymbolicNode[], SYMBOLIC_ONE)]
            with_number(isone, coeff) || (
                terms = distribute(terms, [ExpandedTerm(with_number(constant_node, coeff))])
            )
            for (base, exponent) in dict
                terms = distribute(terms, power_terms(base, exponent))
            end
            terms
        end
        BSImpl.Div(; num, den) => begin
            denominator = expanded_expression(expanded_terms(den))
            ExpandedTerm[
                ExpandedTerm(
                    term.factors, denominator_product(term.denominator, denominator)
                ) for term in expanded_terms(num)
            ]
        end
        BSImpl.Term(; f, args) && if f === (^) && isconst(args[2])
        end => power_terms(args[1], unwrap_const(args[2]))
        _ => [ExpandedTerm(x)]
    end
end

"Append `terms`, each multiplied by the numeric `weight`, onto `destination`."
function append_weighted!(
    destination::Vector{ExpandedTerm}, terms::Vector{ExpandedTerm}, @nospecialize(weight)
)
    if with_number(isone, weight)
        append!(destination, terms)
    else
        factor = with_number(constant_node, weight)
        for term in terms
            push!(
                destination,
                ExpandedTerm(push!(copy(term.factors), factor), term.denominator),
            )
        end
    end
    return destination
end

"The expanded terms of `base^exponent` for a numeric `exponent`."
function power_terms(base::SymbolicNode, @nospecialize(exponent))
    with_number(isone, exponent) && return expanded_terms(base)
    if with_number(is_inverse_power, exponent)
        n = with_number(repetitions, exponent)
        denominator = expanded_expression(power_terms(base, n))
        return [ExpandedTerm(SymbolicNode[], denominator)]
    end
    with_number(is_distributed_power, exponent) ||
        return [ExpandedTerm(base^with_number(constant_node, exponent))]
    n = with_number(repetitions, exponent)
    terms = expanded_terms(base)
    if length(terms) == 1
        term = only(terms)
        isequal(term.denominator, SYMBOLIC_ONE) ||
            return [ExpandedTerm(base^with_number(constant_node, exponent))]
        return [ExpandedTerm(repeat(term.factors, n), SYMBOLIC_ONE)]
    end
    powers = terms
    for _ in 2:n
        powers = distribute(powers, terms)
    end
    return powers
end

"The term as an expression."
function term_expression(term::ExpandedTerm)
    numerator = monomial(term.factors)
    isequal(term.denominator, SYMBOLIC_ONE) && return numerator
    return numerator / term.denominator
end

"The sum of `terms` as an expression, with like terms combined."
function expanded_expression(terms::Vector{ExpandedTerm})
    return polynomial(SymbolicNode[term_expression(term) for term in terms])
end
