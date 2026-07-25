struct BareissFailure <: Exception
    message::String
end

const _BAREISS_VARIABLE = only(DP.@polyvar __questbase_bareiss_variable)

function Base.showerror(io::IO, error::BareissFailure)
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
                augmented[pivot_index, column], augmented[pivot_row, column] =
                    augmented[pivot_row, column], augmented[pivot_index, column]
            end
        end

        old = copy(augmented)
        pivot = old[pivot_index, pivot_index]
        for row in 1:n
            row == pivot_index && continue
            for column in axes(augmented, 2)
                column == pivot_index && continue
                numerator =
                    pivot * old[row, column] -
                    old[row, pivot_index] * old[pivot_index, column]
                augmented[row, column] = if pivot_index == 1
                    numerator
                else
                    _exact_polynomial_division(numerator, previous_pivot)
                end
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

_to_polynomial(expression::Num, symbolic_to_poly, poly_to_symbolic) = _to_polynomial(
    unwrap(expression), symbolic_to_poly, poly_to_symbolic
)
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
        MP.similar_variable(_BAREISS_VARIABLE, name)
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
