"""
Role of a [`HarmonicVariable`](@ref) in the ansatz of its natural variable.

- `CosineQuadrature`: the variable `u` multiplying `cos(ωt)`.
- `SineQuadrature`: the variable `v` multiplying `sin(ωt)`.
- `ZeroHarmonic`: the constant part `a`.
- `LimitCycleFrequency`: the limit-cycle frequency that replaces a gauge-fixed quadrature.
- `Unassigned`: a bare symbol outside any harmonic ansatz.
"""
@enum HarmonicKind begin
    CosineQuadrature
    SineQuadrature
    ZeroHarmonic
    LimitCycleFrequency
    Unassigned
end

"""
$(TYPEDSIGNATURES)

Short label of `kind` used in variable names: `u`, `v`, `a`, `Hopf` or the empty string.
"""
function kind_label(kind::HarmonicKind)
    kind == CosineQuadrature && return "u"
    kind == SineQuadrature && return "v"
    kind == ZeroHarmonic && return "a"
    kind == LimitCycleFrequency && return "Hopf"
    return ""
end

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
    """Role of the variable in the ansatz."""
    kind::HarmonicKind
    """The harmonic being described."""
    ω::Num
    """The natural variable whose harmonic is being described."""
    natural_variable::Num
end

function HarmonicVariable(symbol::Num)
    return HarmonicVariable(symbol, "", Unassigned, Num(1), Num(0))
end

function Base.copy(hv::HarmonicVariable)
    return HarmonicVariable(hv.symbol, hv.name, hv.kind, hv.ω, hv.natural_variable)
end

function Base.show(io::IO, hv::HarmonicVariable)
    if hv.kind == Unassigned
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
function show_ansatz(var::HarmonicVariable)
    var.kind in (CosineQuadrature, SineQuadrature) || return string(var.symbol)
    t = arguments(unwrap(var.natural_variable))
    t = length(t) == 1 ? string(t[1]) : error("more than 1 independent variable")
    trig = var.kind == CosineQuadrature ? "cos" : "sin"
    return string(var.symbol) * "*" * trig * "(" * string(var.ω) * t * ")"
end

Base.display(var::HarmonicVariable) = display(var.name)
Base.display(var::Vector{HarmonicVariable}) = display.(getfield.(var, Symbol("name")))

function coordinate_transform(new_var::Num, ω::Num, t::Num, kind::HarmonicKind)
    kind == CosineQuadrature && return new_var * cos(ω * t)
    kind == SineQuadrature && return new_var * sin(ω * t)
    kind == ZeroHarmonic && return new_var
    throw(ArgumentError("no coordinate transform for a harmonic variable of kind $kind"))
end

"""
$(TYPEDSIGNATURES)

Creates a new harmonic variable and its corresponding transformation rule.
This function takes a natural variable (`nat_var`), a harmonic frequency (`ω`),
an independent variable (`t`), and a [`HarmonicKind`](@ref) (`kind`), and creates a new harmonic variable with
the specified `new_symbol`. It returns a tuple containing the transformation rule and
the new harmonic variable.
"""
function create_harmonic_variable(
    nat_var::Num, ω::Num, t::Num, kind::HarmonicKind; new_symbol::String
)::Tuple{Num,HarmonicVariable}
    new_var = declare_variable(new_symbol, t)
    name = kind_label(kind) * "_{" * var_name(nat_var) * "," * Base.replace(string(ω), "*" => "") * "}"
    rule = coordinate_transform(new_var, ω, t, kind)
    hvar = HarmonicVariable(new_var, name, kind, ω, nat_var)
    return rule, hvar
end

"when HV is used for substitute, substitute its symbol"
function QuestBase.substitute_all(eq::Union{Num,Equation}, rules::Dict{HarmonicVariable})
    return Symbolics.substitute(
        eq, Dict{Num,Any}(zip([hv.symbol for hv in keys(rules)], values(rules)))
    )
end

function QuestBase.substitute_all(var::HarmonicVariable, rules)
    sym, freq = var.symbol, var.ω
    return HarmonicVariable(
        substitute_all(sym, rules),
        var.name,
        var.kind,
        substitute_all(freq, rules),
        var.natural_variable,
    )
end

function QuestBase.substitute_all(vars::Vector{HarmonicVariable}, rules)
    return [substitute_all(var, rules) for var in vars]
end

"Returns the symbols of a `HarmonicVariable`."
Symbolics.get_variables(var::HarmonicVariable)::Num = Num(first(get_variables(var.symbol)))

Base.isequal(v1::HarmonicVariable, v2::HarmonicVariable)::Bool =
    isequal(v1.symbol, v2.symbol)
