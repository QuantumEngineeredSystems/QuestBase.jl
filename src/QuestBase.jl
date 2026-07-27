module QuestBase

using DocStringExtensions
import DynamicPolynomials as DP
import MultivariatePolynomials as MP
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
    unwrap_const,
    unwrap

using Symbolics:
    Symbolics,
    Num,
    wrap,
    get_variables,
    Equation,
    Differential,
    arguments,
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
include("Symbolics/pythagorean.jl")
include("Symbolics/fourier.jl")
include("Symbolics/drop_powers.jl")
include("Symbolics/linear_solve.jl")
include("DifferentialEquation.jl")
include("Variables.jl")
include("HarmonicVariable.jl")
include("HarmonicEquation.jl")

end
