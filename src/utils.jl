function dummy_symbolic_Jacobian(n::Int)::Matrix{Num}
    return Num.(float.(collect(LinearAlgebra.I(n))) .* NaN)
end

flatten(a) = collect(Iterators.flatten(a))

"""Strip Symbolics derivative wrappers to recover the base dependent variable."""
strip_derivative(x::Num) = wrap(strip_derivative(unwrap(x)))
function strip_derivative(x::BasicSymbolic)
    y = x
    while Symbolics.is_derivative(y)
        y = first(arguments(y))
    end
    return y
end
strip_derivative(x) = x

"Show fields of an object."
function show_fields(object)
    for field in fieldnames(typeof(object)) # display every field
        display(string(field))
        display(getfield(object, field))
    end
end

_is_symbolic_like(x) = x isa Num || x isa BasicSymbolic

function _eqtest_symbolic_scalar(a, b)
    diff = Symbolics.simplify(Symbolics.expand(a - b))
    return isequal(diff, 0)
end

function _eqtest_equal(a::Complex, b::Complex)
    _eqtest_equal(real(a), real(b)) && _eqtest_equal(imag(a), imag(b))
end
_eqtest_equal(a::Complex, b) = _eqtest_equal(real(a), b) && _eqtest_equal(imag(a), 0)
_eqtest_equal(a, b::Complex) = _eqtest_equal(a, real(b)) && _eqtest_equal(0, imag(b))

function _eqtest_equal(a::AbstractArray, b::AbstractArray)
    size(a) == size(b) || return false
    return all(i -> _eqtest_equal(a[i], b[i]), eachindex(a, b))
end

function _eqtest_equal(a::Tuple, b::Tuple)
    length(a) == length(b) || return false
    return all(i -> _eqtest_equal(a[i], b[i]), eachindex(a))
end

function _eqtest_equal(a::Equation, b::Equation)
    _eqtest_equal(a.lhs, b.lhs) && _eqtest_equal(a.rhs, b.rhs)
end

function _eqtest_equal(a, b)
    isequal(a, b) && return true
    (_is_symbolic_like(a) || _is_symbolic_like(b)) || return false
    return try
        _eqtest_symbolic_scalar(a, b)
    catch
        false
    end
end

_eqtest_notequal(a, b) = !_eqtest_equal(a, b)

macro eqtest(expr)
    @assert expr.head == :call && expr.args[1] in [:(==), :(!=)]
    qb = QuoteNode(QuestBase)
    return esc(
        if expr.args[1] == :(==)
            :(@test $qb._eqtest_equal($(expr.args[2]), $(expr.args[3])))
        else
            :(@test $qb._eqtest_notequal($(expr.args[2]), $(expr.args[3])))
        end,
    )
end

macro eqsym(expr)
    @assert expr.head == :call && expr.args[1] in [:(==), :(!=)]
    return esc(
        if expr.args[1] == :(==)
            :(isequal($(expr.args[2]), $(expr.args[3])))
        else
            :(!isequal($(expr.args[2]), $(expr.args[3])))
        end,
    )
end

is_identity(A::Matrix{Num}) = (@eqsym A == Matrix{Num}(LinearAlgebra.I, size(A)...))
hasnan(x::Matrix{Num}) = any(my_isnan, Symbolics.value.(x))
my_isnan(x) = isnan(x)
my_isnan(x::BasicSymbolic) = false
