module QuestBase

using DocStringExtensions
using OrderedCollections: OrderedCollections, OrderedDict, OrderedSet
using LinearAlgebra: LinearAlgebra

using SymbolicUtils:
    SymbolicUtils,
    Postwalk,
    BasicSymbolic,
    unwrap,
    isterm,
    ispow,
    isadd,
    isdiv,
    ismul,
    add_with_div,
    issym

using Symbolics:
    Symbolics,
    Num,
    wrap,
    get_variables,
    Equation,
    Differential,
    arguments,
    substitute,
    term,
    expand,
    operation,
    expand_derivatives

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
