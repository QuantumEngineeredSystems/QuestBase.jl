"""
$(TYPEDEF)

Holds a set of algebraic equations governing the harmonics of a system of equations of
motion. The system it was derived from is retained in `source_equations` and is retrievable
with [`source`](@ref); its type is the type parameter `T` and is retrievable with
[`source_type`](@ref).

Typically `T` is a [`DifferentialEquation`](@ref), but a `HarmonicEquation` can also be
derived from a
[`QuantumCumulants.MeanfieldEquations`](https://qojulia.github.io/QuantumCumulants.jl/stable/api/#QuantumCumulants.MeanfieldEquations).

# Fields
$(TYPEDFIELDS)
"""
mutable struct HarmonicEquation{T}
    """A set of equations governing the harmonics."""
    equations::Vector{Equation}
    """A set of variables describing the harmonics."""
    variables::Vector{HarmonicVariable}
    """The parameters of the equation set."""
    parameters::Vector{Num}
    "The Jacobian of the harmonic equations."
    jacobian::Matrix{Num}
    "The system of equations of motion the harmonic equations were derived from."
    source_equations::T

    # use a self-referential constructor with _parameters
    function HarmonicEquation(
        equations::Vector{Equation},
        variables::Vector{HarmonicVariable},
        nat_eq::DifferentialEquation,
    )
        return (
            x=new{DifferentialEquation}(
                equations,
                variables,
                Num[],
                dummy_symbolic_Jacobian(length(variables)),
                nat_eq,
            );
            x.parameters=_parameters(x);
            x
        )
    end
    function HarmonicEquation(
        equations::Vector{Equation},
        variables::Vector{HarmonicVariable},
        parameters::Vector{Num},
        source_equations::DifferentialEquation,
    )
        return new{DifferentialEquation}(
            equations,
            variables,
            parameters,
            dummy_symbolic_Jacobian(length(variables)),
            source_equations,
        )
    end
    function HarmonicEquation(
        equations::Vector{Equation},
        variables::Vector{HarmonicVariable},
        parameters::Vector{Num},
        jacobian::Matrix{Num},
        source_equations::T,
    ) where {T}
        return new{T}(equations, variables, parameters, jacobian, source_equations)
    end
end

"""
$(TYPEDSIGNATURES)

Return the system of equations of motion `eom` was derived from, e.g. the
[`DifferentialEquation`](@ref) the harmonic ansatz was applied to.

See also [`source_type`](@ref).
"""
source(eom::HarmonicEquation) = eom.source_equations

"""
$(TYPEDSIGNATURES)

Return the type of the system of equations of motion `eom` was derived from. Use this to
dispatch on the origin of a `HarmonicEquation` without materialising the source system.

See also [`source`](@ref).
"""
source_type(eom::HarmonicEquation{T}) where {T} = T

"Get the parameters (not time nor variables) of a HarmonicEquation"
function _parameters(eom::HarmonicEquation)
    all_symbols = Num[]
    for eq in eom.equations
        for v in collect(Symbolics.get_variables(eq.lhs))
            is_derivative(v) || push!(all_symbols, Num(v))
        end
        for v in collect(Symbolics.get_variables(eq.rhs))
            is_derivative(v) || push!(all_symbols, Num(v))
        end
    end
    # subtract the set of independent variables (i.e., time) from all free symbols
    return unique(setdiff(all_symbols, get_variables(eom), get_independent_variables(eom)))
end

"""
$(TYPEDSIGNATURES)
Get the internal symbols of the independent variables of `eom`.
"""
function Symbolics.get_variables(eom::HarmonicEquation)::Vector{Num}
    return get_variables.(eom.variables)
end

function Base.show(io::IO, eom::HarmonicEquation)
    println(io, "A set of ", length(eom.equations), " harmonic equations")
    println(io, "Variables: ", join(string.(get_variables(eom)), ", "))
    println(io, "Parameters: ", join(string.(eom.parameters), ", "))
    println(io, "\nHarmonic ansatz: ", _show_ansatz(eom))
    println(io, "\nHarmonic equations:")
    return [println(io, "\n", eq) for eq in eom.equations]
end

"""Gives the full harmonic ansatz used to construct `eom`."""
function _show_ansatz(eom::HarmonicEquation)
    output = ""
    vars = unique(getfield.(eom.variables, :natural_variable))
    for nat_var in vars
        # the Hopf variable (limit cycle frequency) does not contribute a term
        harm_vars = filter(
            x -> isequal(nat_var, x.natural_variable) && x.type !== "Hopf", eom.variables
        )
        ansatz = join([_show_ansatz(var) for var in harm_vars], " + ")
        output *= "\n" * string(nat_var) * " = " * ansatz
    end
    return output
end

Base.show(eom::HarmonicEquation) = show_fields(eom)

"Apply `rules` to both `equations` and `variables` field of `eom`"
function substitute_all(eom::HarmonicEquation, rules::Union{Dict,Pair})::HarmonicEquation
    new_eom = deepcopy(eom)
    new_eom.equations = expand_derivatives.(substitute_all(eom.equations, rules))
    return new_eom
end

#   Drop powers of `var` of degree >= `deg` from the equation set in `eom`.
function drop_powers(eom::HarmonicEquation, terms::Vector{Num}, deg::Int)
    new_eom = deepcopy(eom)
    new_eom.equations = drop_powers(eom.equations, terms, deg)
    return new_eom
end

"""
$(TYPEDSIGNATURES)
Return the independent variables (typically time) of `eom`.
"""
function get_independent_variables(eom::HarmonicEquation)::Vector{Num}
    dynamic_vars = flatten(getfield.(eom.variables, Symbol("symbol")))
    return flatten(unique([arguments(unwrap(var)) for var in dynamic_vars]))
end

_remove_brackets(var::Num) = declare_variable(var_name(var))
_remove_brackets(vars::Vector{Num}) = _remove_brackets.(vars)
_remove_brackets(hv::HarmonicVariable) = _remove_brackets(hv.symbol)

"Returns the equation system in `eom`, dropping all argument brackets (i.e., u(T) becomes u)."
function _remove_brackets(eom::HarmonicEquation)
    variable_rules = [var => _remove_brackets(var) for var in get_variables(eom)]
    equations_lhs = Num.(getfield.(eom.equations, :lhs) - getfield.(eom.equations, :rhs))
    return substitute_all(equations_lhs, variable_rules)
end

function _remove_brackets(eqs::Vector{Num}, vars::Vector{Num})
    vars_ = _remove_brackets.(vars)
    variable_rules = Dict(zip(vars, vars_))
    return substitute_all(eqs, variable_rules), vars_
end

"""
$(TYPEDSIGNATURES)
Rearrange `eom` to the standard form, such that the derivatives of the variables are on one side.
"""
function rearrange_standard(eom::HarmonicEquation)
    tvar = get_independent_variables(eom)[1]
    dvars = d(get_variables(eom), tvar)
    return rearrange(eom, dvars)
end

"""
$(TYPEDSIGNATURES)
Rearrange `eom` to the standard form, such that the derivatives of the variables are on one side.
"""
function rearrange_standard!(eom::HarmonicEquation)
    tmp_eom = rearrange_standard(eom::HarmonicEquation)
    eom.equations = tmp_eom.equations
    return eom
end

"Rearrange an equation system such that the field equations is equal to the vector specified in new_lhs"
function rearrange!(eom::HarmonicEquation, new_rhs::Vector{Num})
    soln = try
        fraction_free_linear_solve(eom.equations, new_rhs)
    catch error
        error isa BareissFailure || rethrow()
        fallback = Symbolics.symbolic_linear_solve(
            eom.equations, new_rhs; simplify=false, check=true
        )
        Symbolics.simplify_fractions.(Num.(fallback))
    end
    # SymbolicUtils 4 leaves `cos(x)^2 + sin(x)^2` standing, and these solutions carry the
    # coefficient determinant in their denominator, so without collapsing it the expressions
    # grow every time a system is rearranged.
    soln = collapse_pythagorean.(soln)
    soln = reduce_denominator.(soln)
    eom.equations = soln .~ new_rhs
    return nothing
end

"Rearrange an equation system such that the field equations is equal to the vector specified in new_lhs"
function rearrange(eom::HarmonicEquation, new_rhs::Vector{Num})
    new_eom = deepcopy(eom)
    rearrange!(new_eom, new_rhs)
    return new_eom
end

"""
$(TYPEDSIGNATURES)
Check if `eom` is rearranged to the standard form, such that the derivatives of the variables are on one side.
"""
function is_rearranged(eom::HarmonicEquation)
    tvar = get_independent_variables(eom)[1]
    dvar = d(get_variables(eom), tvar)
    lhs = getfield.(eom.equations, :lhs)
    rhs = getfield.(eom.equations, :rhs)

    HB_bool = isequal(rhs, dvar)
    hopf_bool = in("Hopf", getfield.(eom.variables, :type))
    # no derivative anywhere on the left-hand side, asked of the tree rather than `string(lhs)`
    MF_bool = !any(occurs_in(dv, l) for dv in dvar, l in lhs)

    # Hopf-containing equations or MF equation are arranged by construstion
    return HB_bool || hopf_bool || MF_bool
end

"""
$(TYPEDSIGNATURES)

(Re)-declare the variables of `eom` in module scope.
"""
function declare_variables(eom::HarmonicEquation)
    vars_orig = get_variables(eom)
    return declare_variable.(var_name.(vars_orig))
end
# TODO should this be the variable with independent_variable or without independent_variable?
