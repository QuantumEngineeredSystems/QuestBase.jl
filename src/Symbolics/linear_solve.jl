struct BareissFailure <: Exception
    message::String
end

"""
Raised when Bareiss elimination is abandoned because its intermediate polynomials grew past
what is workable, rather than because the system is unsuitable for it.

This is deliberately *not* a [`BareissFailure`](@ref). A `BareissFailure` means the system is
not polynomial in the unknowns, and `rearrange!` answers it by falling back to
`Symbolics.symbolic_linear_solve`. Size is different: the LU fallback grows nested fractions
on exactly the systems that overwhelm Bareiss, so it is no refuge. This propagates to the
caller instead, letting one that has a cheaper route (`add_jacobian!` evaluates the Jacobian
implicitly) take it.
"""
struct BareissTooLarge <: Exception
    message::String
end

const _BAREISS_VARIABLE = only(DP.@polyvar __questbase_bareiss_variable)

"""
Largest number of monomials an intermediate Bareiss entry may reach before the solve gives
up with a [`BareissFailure`](@ref).

Bareiss is fraction-free, which is what lets it sidestep the nested-fraction growth of
`Symbolics.symbolic_linear_solve`, but the entries it carries are exact polynomials whose
degree grows with every pivot. On a dense system in many symbolic atoms, such as the
harmonic equations of a van der Pol oscillator expanded in three harmonics, that growth is
combinatorial and the solve never finishes even though nothing has gone wrong. Callers that
have a cheaper option, such as evaluating a Jacobian implicitly, can catch the failure and
take it.

The systems this solve handles well stay far below the budget: 10 monomials for a Duffing
oscillator, 28 for a parametron, 344 for a van der Pol with two harmonics. The budget has to
be tight as well as safe, because the polynomial arithmetic is already slow at entry sizes
well under it, so a generous budget takes minutes to trip.
"""
const BAREISS_TERM_BUDGET = Ref(2_000)

"""
Ceiling on the monomial count a single polynomial product may be predicted to reach.

The size checks on finished entries are not enough on their own: one multiplication of two
wide polynomials can exhaust memory before its result can be measured. The number of terms
in a product is at most the product of the factors' term counts, so that bound is checked
before multiplying. It is deliberately far looser than [`BAREISS_TERM_BUDGET`](@ref),
because the bound is pessimistic and systems that solve fine reach it: a van der Pol with
two harmonics multiplies 344-term entries, predicting over a hundred thousand terms for a
result that is much smaller.
"""
const BAREISS_PRODUCT_BUDGET = Ref(4_000_000)

"Number of monomials in a Bareiss entry; constants count as one."
_bareiss_size(entry) = entry isa Number ? 1 : length(MP.terms(entry))

"Multiply two Bareiss entries, refusing products predicted to be unmanageably wide."
function _bareiss_product(left, right, pivot_index, n)
    predicted = _bareiss_size(left) * _bareiss_size(right)
    predicted > BAREISS_PRODUCT_BUDGET[] && throw(
        BareissTooLarge(
            "Bareiss product predicted at $(predicted) terms at pivot $(pivot_index) of " *
            "$(n), over the $(BAREISS_PRODUCT_BUDGET[]) term product budget",
        ),
    )
    return left * right
end

function Base.showerror(io::IO, error::BareissFailure)
    return print(io, error.message)
end

function Base.showerror(io::IO, error::BareissTooLarge)
    return print(io, error.message)
end

"""
    fraction_free_linear_solve(equations, variables)

Solve a square symbolic linear system using pivoted Bareiss–Jordan elimination.

The coefficients are converted together into SymbolicUtils' sparse polynomial
representation. Bareiss elimination then performs only exact polynomial divisions
and produces one determinant denominator per solution, instead of the nested
fractions produced by symbolic LU.
"""
function fraction_free_linear_solve(equations, variables)
    coefficients, offsets, islinear = Symbolics.linear_expansion(equations, variables)
    islinear || throw(ArgumentError("the equation system is not linear in the variables"))

    rows, columns = size(coefficients)
    rows == columns == length(offsets) ||
        throw(DimensionMismatch("fraction-free linear solve requires a square system"))

    return _bareiss_jordan_solve(coefficients, -offsets)
end

function _bareiss_jordan_solve(coefficients, right_hand_side)
    n = size(coefficients, 1)
    symbolic_type = eltype(coefficients)
    poly_variable_type = typeof(_BAREISS_VARIABLE)
    poly_to_symbolic = Dict{poly_variable_type,symbolic_type}()
    symbolic_to_poly = Dict{symbolic_type,poly_variable_type}()
    augmented = Matrix{Any}(undef, n, n + 1)

    for column in 1:n, row in 1:n
        augmented[row, column] = _to_polynomial(
            coefficients[row, column], symbolic_to_poly, poly_to_symbolic
        )
    end
    for row in 1:n
        augmented[row, n + 1] = _to_polynomial(
            right_hand_side[row], symbolic_to_poly, poly_to_symbolic
        )
    end

    previous_pivot = 1
    for pivot_index in 1:n
        pivot_row = _find_bareiss_pivot(augmented, pivot_index, n)
        isnothing(pivot_row) && throw(LinearAlgebra.SingularException(pivot_index))
        if pivot_row != pivot_index
            for column in axes(augmented, 2)
                swapped = augmented[pivot_index, column]
                augmented[pivot_index, column] = augmented[pivot_row, column]
                augmented[pivot_row, column] = swapped
            end
        end

        old = copy(augmented)
        pivot = old[pivot_index, pivot_index]
        for row in 1:n
            row == pivot_index && continue
            for column in axes(augmented, 2)
                column == pivot_index && continue
                numerator =
                    _bareiss_product(pivot, old[row, column], pivot_index, n) -
                    _bareiss_product(
                        old[row, pivot_index], old[pivot_index, column], pivot_index, n
                    )
                entry = if pivot_index == 1
                    numerator
                else
                    _exact_polynomial_division(numerator, previous_pivot)
                end
                terms = _bareiss_size(entry)
                terms > BAREISS_TERM_BUDGET[] && throw(
                    BareissTooLarge(
                        "Bareiss entry reached $(terms) terms at pivot $(pivot_index) of " *
                        "$(n), over the $(BAREISS_TERM_BUDGET[]) term budget",
                    ),
                )
                augmented[row, column] = entry
            end
            augmented[row, pivot_index] = 0
        end
        previous_pivot = pivot
    end

    solution = Vector{Num}(undef, n)
    for row in 1:n
        numerator_poly = augmented[row, n + 1]
        denominator_poly = augmented[row, row]
        leading_coefficient = if denominator_poly isa Number
            denominator_poly
        else
            MP.leading_coefficient(denominator_poly)
        end
        if leading_coefficient isa Real && leading_coefficient < 0
            numerator_poly = -numerator_poly
            denominator_poly = -denominator_poly
        end
        numerator = _from_polynomial(numerator_poly, poly_to_symbolic)
        denominator = _from_polynomial(denominator_poly, poly_to_symbolic)
        solution[row] = numerator / denominator
    end
    return solution
end

function _to_polynomial(expression::Num, symbolic_to_poly, poly_to_symbolic)
    return _to_polynomial(unwrap(expression), symbolic_to_poly, poly_to_symbolic)
end
_to_polynomial(expression::Number, _, _) = expression
function _to_polynomial(expression::BasicSymbolic, symbolic_to_poly, poly_to_symbolic)
    if SymbolicUtils.isconst(expression)
        return unwrap_const(expression)
    elseif isadd(expression)
        result = 0
        for argument in arguments(expression)
            result += _to_polynomial(argument, symbolic_to_poly, poly_to_symbolic)
        end
        return result
    elseif ismul(expression)
        result = 1
        for argument in arguments(expression)
            result *= _to_polynomial(argument, symbolic_to_poly, poly_to_symbolic)
        end
        return result
    elseif ispow(expression)
        base, exponent = arguments(expression)
        if SymbolicUtils.isconst(exponent)
            value = unwrap_const(exponent)
            if value isa Integer && value >= 0
                return _to_polynomial(base, symbolic_to_poly, poly_to_symbolic)^value
            end
        end
    end

    variable = get!(symbolic_to_poly, expression) do
        name = Symbol(:questbase_atom_, hash(expression))
        return MP.similar_variable(_BAREISS_VARIABLE, name)
    end
    get!(poly_to_symbolic, variable, expression)
    return variable
end

_from_polynomial(value::Number, _) = Num(value)
function _from_polynomial(variable::typeof(_BAREISS_VARIABLE), poly_to_symbolic)
    return Num(poly_to_symbolic[variable])
end
function _from_polynomial(polynomial, poly_to_symbolic)
    variables = MP.variables(polynomial)
    result = Num(0)
    for term in MP.terms(polynomial)
        expression = Num(MP.coefficient(term))
        exponents = MP.exponents(MP.monomial(term))
        for (variable, exponent) in zip(variables, exponents)
            iszero(exponent) && continue
            expression *= Num(poly_to_symbolic[variable])^exponent
        end
        result += expression
    end
    return result
end

function _find_bareiss_pivot(augmented, pivot_index, n)
    best_row = nothing
    best_size = typemax(Int)
    for row in pivot_index:n
        candidate = augmented[row, pivot_index]
        iszero(candidate) && continue
        candidate_size = candidate isa Number ? 1 : length(MP.terms(candidate))
        if candidate_size < best_size
            best_row = row
            best_size = candidate_size
        end
    end
    return best_row
end

function _exact_polynomial_division(numerator, denominator)
    quotient = try
        MP.div_multiple(numerator, denominator)
    catch error
        throw(BareissFailure("Bareiss exact division failed: $(sprint(showerror, error))"))
    end
    iszero(quotient * denominator - numerator) ||
        throw(BareissFailure("Bareiss division produced a nonzero remainder"))
    return quotient
end
