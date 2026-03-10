module QuestBase

using DocStringExtensions
using OrderedCollections: OrderedCollections, OrderedDict, OrderedSet
using LinearAlgebra: LinearAlgebra

using SymbolicUtils:
    SymbolicUtils,
    Postwalk,
    BasicSymbolic,
    isterm,
    ispow,
    isadd,
    isdiv,
    ismul,
    issym,
    add_with_div,
    is_literal_number,
    unwrap_const

using Symbolics:
    Symbolics,
    Num,
    unwrap,
    wrap,
    get_variables,
    Equation,
    Differential,
    arguments,
    substitute,
    term,
    expand,
    operation,
    expand_derivatives,
    tosymbol,
    is_derivative

include("utils.jl")
include("Symbolics/Symbolics_utils.jl")
include("Symbolics/exponentials.jl")
include("Symbolics/fourier.jl")
include("Symbolics/drop_powers.jl")
include("DifferentialEquation.jl")
include("Variables.jl")
include("HarmonicVariable.jl")
include("HarmonicEquation.jl")

end
