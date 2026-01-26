"The derivative of f w.r.t. x of degree deg"
function d(f::Num, x::Num, deg=1)::Num
    return isequal(deg, 0) ? f : (Differential(x)^deg)(f)
end
d(funcs::Vector{Num}, x::Num, deg=1) = Num[d(f, x, deg) for f in funcs]

"Declare a variable in the the currently active Module namespace"
function declare_variable(name::String)
    var_sym = Symbol(name)
    @eval($(var_sym) = first(Symbolics.@variables $var_sym))
    return eval(var_sym)
end

declare_variable(x::Num) = declare_variable(string(x))

"Declare a variable that is a function of another variable in the Module namespace"
function declare_variable(name::String, independent_variable::Num)
    # independent_variable = declare_variable(independent_variable) convert string into Num
    var_sym = Symbol(name)
    new_var = Symbolics.@variables $var_sym(independent_variable)
    @eval($(var_sym) = first($new_var)) # store the variable under "name" in this namespace
    return eval(var_sym)
end

"Return the name of a variable (excluding independent variables)"
function var_name(x::Num)::String
    var = string(x)
    var = replace(var, r"\(.*\)$" => "")
    return String(replace(var, r"\\mathtt\{([^}]*)\}" => s"\1"))
    # ^ remove "\\mathtt{}" from the variable name coming from Symbolics
    # since Symbolics v6.14.1 (Symbolics#1305)
end
var_name(x::SymbolicUtils.Sym) = String(x.name)
