"The derivative of f w.r.t. x of degree deg"
function d(f::Num, x::Num, deg=1)::Num
    return isequal(deg, 0) ? f : (Differential(x)^deg)(f)
end
d(funcs::Vector{Num}, x::Num, deg=1) = Num[d(f, x, deg) for f in funcs]

"""
The symtype `Symbolics.@variables f(t)` gives its variable. Spelled out in full because
`Symbolics.variable` needs a concrete type and `FnType{Tuple,Real}` is still a `UnionAll`.
"""
const _FnTypeReal = SymbolicUtils.FnType{Tuple,Real,Nothing}

"""
$(TYPEDSIGNATURES)

Declare a symbolic variable named `name`.

`Symbolics.@variables` needs the name as a literal, so names only known at runtime (`u1`,
`v2`, the bracket-free copies [`remove_brackets`](@ref) makes) have to go through
`Symbolics.variable` instead.

This used to build the variable with `@eval` and bind the result inside this module. Nothing
ever read that binding, since it lands in `QuestBase` rather than in the caller's namespace,
and evaluating into a closed module makes the package impossible to precompile with a
workload: `@compile_workload` runs the very code that creates these variables.
"""
declare_variable(name::String) = Symbolics.variable(Symbol(name))

declare_variable(x::Num) = declare_variable(string(x))

"""
$(TYPEDSIGNATURES)

Declare a variable named `name` that is a function of `independent_variable`.
"""
function declare_variable(name::String, independent_variable::Num)
    function_symbol = SymbolicUtils.Sym{NUM_VARTYPE}(Symbol(name); type=_FnTypeReal)
    return Num(function_symbol(unwrap(independent_variable)))
end

"Return the name of a variable (excluding independent variables)"
var_name(x::Num) = var_name(unwrap(x))
function var_name(x::SymbolicNode)
    return @match x begin
        BSImpl.Sym(; name) => String(name)
        BSImpl.Term(; f) && if f isa SymbolicNode
        end => var_name(f)
        _ => throw(ArgumentError("$x is not a variable."))
    end
end
