"""
$(TYPEDEF)

Holds a variable stored under `symbol` describing the harmonic `ω` of `natural_variable`.

# Fields
$(TYPEDFIELDS)
"""
mutable struct HarmonicVariable
    """Symbol of the variable in the HarmonicBalance namespace."""
    symbol::Num
    """Human-readable labels of the variable, used for plotting."""
    name::String
    """Type of the variable (u or v for quadratures, a for a constant, Hopf for Hopf etc.)"""
    type::String
    """The harmonic being described."""
    ω::Num
    """The natural variable whose harmonic is being described."""
    natural_variable::Num
end

function HarmonicVariable(symbol::Num)
    return HarmonicVariable(symbol, "", "", Num(1), Num(0))
end

function Base.show(io::IO, hv::HarmonicVariable)
    if isempty(hv.type)
        s = "Harmonic variable " * string.(hv.symbol)
    else
        s =
            "Harmonic variable " *
            string.(hv.symbol) *
            " for harmonic " *
            string(hv.ω) *
            " of " *
            string(hv.natural_variable)
    end
    return println(io, s)
end

"""Gives the relation between `var` and the underlying natural variable."""
function _show_ansatz(var::HarmonicVariable)
    if isempty(var.type)
        return string(var.symbol)
    end
    t = arguments(var.natural_variable.val)
    t = length(t) == 1 ? string(t[1]) : error("more than 1 independent variable")
    ω = string(var.ω)
    terms = Dict("u" => "*cos(" * ω * t * ")", "v" => "*sin(" * ω * t * ")", "a" => "")
    return string(string(var.symbol) * terms[var.type])
end

# pretty-printing
Base.display(var::HarmonicVariable) = display(var.name)
Base.display(var::Vector{HarmonicVariable}) = display.(getfield.(var, Symbol("name")))

function _coordinate_transform(new_var, ω, t, type)::Num
    coords = Dict([
        "u" => new_var * cos(ω * t), "v" => new_var * sin(ω * t), "a" => new_var
    ])
    return coords[type]
end

"""
$(TYPEDSIGNATURES)

Creates a new harmonic variable and its corresponding transformation rule.
This function takes a natural variable (`nat_var`), a harmonic frequency (`ω`),
an independent variable (`t`), and a type (`type`), and creates a new harmonic variable with
the specified `new_symbol`. It returns a tuple containing the transformation rule and
the new harmonic variable.
"""
function _create_harmonic_variable(
    nat_var::Num, ω::Num, t::Num, type::String; new_symbol::String
)::Tuple{Num,HarmonicVariable}
    new_var = declare_variable(new_symbol, t) # this holds the internal symbol
    name = type * "_{" * var_name(nat_var) * "," * Base.replace(string(ω), "*" => "") * "}"
    # contribution of this harmonic variable to the natural variable
    rule = _coordinate_transform(new_var, ω, t, type)
    hvar = HarmonicVariable(new_var, name, type, ω, nat_var)
    return rule, hvar
end

###
# Functions for variable substutions and manipulation of HarmonicVariable
###

"when HV is used for substitute, substitute its symbol"
function QuestBase.substitute_all(eq::Union{Num,Equation}, rules::Dict{HarmonicVariable})
    return Symbolics.substitute(
        eq, Dict(zip(getfield.(keys(rules), :symbol), values(rules)))
    )
end

function QuestBase.substitute_all(var::HarmonicVariable, rules)
    sym, freq = var.symbol, var.ω
    return HarmonicVariable(
        substitute_all(sym, rules),
        var.name,
        var.type,
        substitute_all(freq, rules),
        var.natural_variable,
    )
end

function QuestBase.substitute_all(vars::Vector{HarmonicVariable}, rules)
    return [substitute_all(var, rules) for var in vars]
end

"Returns the symbols of a `HarmonicVariable`."
Symbolics.get_variables(var::HarmonicVariable)::Num = Num(first(get_variables(var.symbol)))

Base.isequal(v1::HarmonicVariable, v2::HarmonicVariable)::Bool = isequal(
    v1.symbol, v2.symbol
)
