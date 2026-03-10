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
    issym,
    add_with_div,
    is_literal_number,
    unwrap_const,
    unwrap,
    AddMulVariant

using Moshi.Match: @match
const BSImpl = SymbolicUtils.BasicSymbolicImpl

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
    expand_derivatives,
    tosymbol,
    is_derivative,
    sorted_arguments

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
