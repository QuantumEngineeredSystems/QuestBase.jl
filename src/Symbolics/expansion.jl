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

"""
    expanded_terms(x)

The terms of `x` with every product and positive integer power of a sum distributed, and
every fraction split over the terms of its numerator. Denominators are expanded on their
own and kept whole. Arguments of function calls are left as they stand.
"""
function expanded_terms(x::SymbolicNode)
    return @match x begin
        BSImpl.AddMul(; variant) => if variant === AddMulVariant.ADD
            terms = ExpandedTerm[]
            for argument in arguments(x)
                append!(terms, expanded_terms(argument))
            end
            terms
        else
            terms = [ExpandedTerm(SymbolicNode[], SYMBOLIC_ONE)]
            for argument in arguments(x)
                terms = distribute(terms, expanded_terms(argument))
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
        BSImpl.Term(; f, args) &&
        if f === (^) && isconst(args[2]) && with_const(is_distributed_power, args[2])
        end => begin
            base = expanded_terms(args[1])
            length(base) == 1 && return [ExpandedTerm(x)]
            terms = base
            for _ in 2:with_exact(Int, args[2])
                terms = distribute(terms, base)
            end
            terms
        end
        _ => [ExpandedTerm(x)]
    end
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
