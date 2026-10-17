"""
$(TYPEDEF)

Holds differential equation(s) of motion and a set of harmonics to expand each variable.
This is the primary input for `HarmonicBalance.jl`. After inputting the equations,
the harmonics ansatz needs to be specified using `add_harmonic!`.

# Fields
$(TYPEDFIELDS)

## Example
```julia-repl
julia> @variables t, x(t), y(t), ω0, ω, F, k;

# equivalent ways to enter the simple harmonic oscillator
julia> DifferentialEquation(d(x,t,2) + ω0^2 * x - F * cos(ω*t), x);
julia> DifferentialEquation(d(x,t,2) + ω0^2 * x ~ F * cos(ω*t), x);

# two coupled oscillators, one of them driven
julia> DifferentialEquation(
    [d(x,t,2) + ω0^2 * x - k*y, d(y,t,2) + ω0^2 * y - k*x] .~ [F * cos(ω*t), 0], [x,y]
);
```
"""
mutable struct DifferentialEquation
    """Assigns to each variable an equation of motion."""
    equations::OrderedDict{Num,Equation}
    """Assigns to each variable a set of harmonics."""
    harmonics::OrderedDict{Num,OrderedSet{Num}}

    function DifferentialEquation(eqs)
        return new(eqs, OrderedDict(var => OrderedSet() for var in keys(eqs)))
    end

    function DifferentialEquation()
        return new(OrderedDict{Num,Equation}(), OrderedDict{Num,OrderedSet{Num}}())
    end

    function DifferentialEquation(eqs::Vector{Equation}, vars::Vector{Num})
        if length(eqs) != length(vars)
            throw(ArgumentError("The number of equations and variables are not the same."))
        end
        return DifferentialEquation(OrderedDict(zip(vars, eqs)))
    end

    function DifferentialEquation(exprs::Vector{Num}, vars::Vector{Num})
        return DifferentialEquation(exprs .~ Int(0), vars)
    end

    function DifferentialEquation(eq::Equation, var::Num)
        typerhs = typeof(eq.rhs)
        typelhs = typeof(eq.lhs)
        if eq.rhs isa AbstractVector || eq.lhs isa AbstractVector
            throw(
                ArgumentError(
                    "The equation is of the form $(typerhs)~$(typelhs) is not supported.
                    Commenly one forgot to broadcast the equation symbol `~`."
                ),
            )
        end
        return DifferentialEquation([eq], [var])
    end
    function DifferentialEquation(eq::Equation, vars::Vector{Num})
        typerhs = typeof(eq.rhs)
        typelhs = typeof(eq.lhs)
        return throw(
            ArgumentError(
                "The variables are of type $(typeof(vars)). Whereas you gave one equation of
                type $(typerhs)~$(typelhs). Commenly one forgot to broadcast the equation symbol `~`.",
            ),
        )
    end
    DifferentialEquation(lhs::Num, var::Num) = DifferentialEquation([lhs ~ Int(0)], [var])
end

function Base.copy(diff_eq::DifferentialEquation)
    new_eq = DifferentialEquation()
    new_eq.equations = copy(diff_eq.equations)
    new_eq.harmonics = OrderedDict{Num,OrderedSet{Num}}(
        var => copy(harmonics) for (var, harmonics) in diff_eq.harmonics
    )
    return new_eq
end

"show method of the type `DifferentialEquation`"
function Base.show(io::IO, diff_eq::DifferentialEquation)
    println(io, "System of ", length(keys(diff_eq.equations)), " differential equations")
    println(io, "Variables:       ", join(keys(diff_eq.equations), ", "))
    print(io, "Harmonic ansatz: ")
    for var in keys(diff_eq.harmonics)
        print(io, string(var), " => ", join(string.(diff_eq.harmonics[var]), ", "))
        print(io, ";   ")
    end
    println(io, "\n")
    return [println(io, eq) for eq in values(diff_eq.equations)]
end

"""
$(TYPEDSIGNATURES)

Return the dependent variables of `diff_eom`.
"""
Symbolics.get_variables(diff_eom::DifferentialEquation)::Vector{Num} =
    collect(keys(diff_eom.equations))

"""
$(TYPEDSIGNATURES)

Check if all equations in `diff_eom` are harmonic with respect to `t`. The function takes a
differential equation system `diff_eom` and a variable `t`, and returns `true` if all equations
are harmonic with respect to `t`, otherwise it returns `false`.
"""
QuestBase.is_harmonic(diff_eom::DifferentialEquation, t::Num)::Bool =
    all([is_harmonic(eq, t) for eq in values(diff_eom.equations)])

"""
$(TYPEDSIGNATURES)

Return the terms of `diff_eom` which are not polynomial in its variables, see
[`nonpolynomial_terms`](@ref). An empty result means every equation is a polynomial in the
variables and their derivatives, which averaging requires.
"""
function nonpolynomial_terms(diff_eom::DifferentialEquation)
    vars = get_variables(diff_eom)
    return unique(
        reduce(
            vcat,
            [nonpolynomial_terms(eq, vars) for eq in values(diff_eom.equations)];
            init=BasicSymbolic[],
        ),
    )
end

"""
$(TYPEDSIGNATURES)

Return the independent dependent variables of `diff_eom`.
"""
function get_independent_variables(diff_eom::DifferentialEquation)
    return Num.(flatten(unique([arguments(unwrap(x)) for x in keys(diff_eom.equations)])))
end

"""
$(TYPEDSIGNATURES)

Add the harmonic `ω` to the harmonic ansatz used to expand the variable `var` in `diff_eom`.

## Example

# define the simple harmonic oscillator and specify that x(t) oscillates with frequency ω
```julia-repl
julia> @variables t, x(t), y(t), ω0, ω, F, k;
julia> diff_eq = DifferentialEquation(d(x,t,2) + ω0^2 * x ~ F * cos(ω*t), x);
julia> add_harmonic!(diff_eq, x, ω) # expand x using ω

System of 1 differential equations
Variables:       x(t)
Harmonic ansatz: x(t) => ω;

(ω0^2)*x(t) + Differential(t)(Differential(t)(x(t))) ~ F*cos(t*ω)
```
"""
function add_harmonic!(diff_eom::DifferentialEquation, var::Num, ω)
    push!.(Ref(diff_eom.harmonics[var]), ω)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the equations of `eom`.
"""
get_equations(eom::DifferentialEquation) = collect(values(eom.equations))

"""
$(TYPEDSIGNATURES)

Checks if the differential equations in `eom` are arranged in standard form, where the highest
derivative of each variable appears isolated on the left-hand side. The default degree is 2,
corresponding to second-order differential equations.
"""
function is_rearranged_standard(eom::DifferentialEquation, degree=2)
    tvar = get_independent_variables(eom)[1]
    D = Differential(tvar)^degree
    return isequal(getfield.(values(eom.equations), :lhs), D.(get_variables(eom)))
end

"""
$(TYPEDSIGNATURES)

Rearrange the equations in `eom` such that the expressions in `new_lhs` appear isolated on the
left-hand sides. Modifies `eom` in place and returns it.
"""
function rearrange!(eom::DifferentialEquation, new_lhs::Vector{Num})
    soln = solve_linear_system(get_equations(eom), new_lhs)
    eom.equations = OrderedDict(zip(collect(keys(eom.equations)), new_lhs .~ soln))
    return eom
end

"""
$(TYPEDSIGNATURES)

Non-mutating version of [`rearrange!`](@ref).
"""
function rearrange(eom::DifferentialEquation, new_lhs::Vector{Num})
    return rearrange!(copy(eom), new_lhs)
end

"""
$(TYPEDSIGNATURES)

Rearrange `eom` to standard form, where the highest derivative of each variable (specified by
`degree`, default 2) appears isolated on the left-hand side. Modifies `eom` in place and
returns it.
"""
function rearrange_standard!(eom::DifferentialEquation, degree=2)
    tvar = get_independent_variables(eom)[1]
    D = Differential(tvar)^degree
    return rearrange!(eom, D.(get_variables(eom)))
end

"""
$(TYPEDSIGNATURES)

Non-mutating version of [`rearrange_standard!`](@ref).
"""
function rearrange_standard(eom::DifferentialEquation, degree=2)
    return rearrange_standard!(copy(eom), degree)
end

function get_variables_nums(vars::Vector{Num})
    return unique(flatten([Num.(get_variables(x)) for x in vars]))
end
