function dummy_symbolic_Jacobian(n::Int)::Matrix{Num}
    return Num.(float.(collect(LinearAlgebra.I(n))) .* NaN)
end

flatten(a) = collect(Iterators.flatten(a))

"Show fields of an object."
function show_fields(object)
    for field in fieldnames(typeof(object)) # display every field
        display(string(field))
        display(getfield(object, field))
    end
end

macro eqtest(expr)
    @assert expr.head == :call && expr.args[1] in [:(==), :(!=)]
    return esc(
        if expr.args[1] == :(==)
            :(@test isequal($(expr.args[2]), $(expr.args[3])))
        else
            :(@test !isequal($(expr.args[2]), $(expr.args[3])))
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
hasnan(x::Matrix{Num}) = any(my_isnan, unwrap.(x))
my_isnan(x) = isnan(x)
my_isnan(x::BasicSymbolic) = false
