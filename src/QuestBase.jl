module QuestBase

using DocStringExtensions
using Nemo: Nemo, QQ, QQFieldElem, QQMPolyRing, QQMPolyRingElem
using OrderedCollections: OrderedCollections, OrderedDict, OrderedSet
using LinearAlgebra: LinearAlgebra
using Moshi.Match: @match

using SymbolicUtils:
    SymbolicUtils,
    AddMulVariant,
    BSImpl,
    Postwalk,
    BasicSymbolic,
    isadd,
    ismul,
    issym,
    unwrap_const,
    unwrap,
    vartype,
    add_worker,
    mul_worker,
    Const,
    isconst,
    iscall

using Symbolics:
    Symbolics,
    Num,
    wrap,
    get_variables,
    Equation,
    Differential,
    arguments,
    term,
    operation,
    expand_derivatives,
    is_derivative

include("utils.jl")
include("Symbolics/Symbolics_utils.jl")
include("Symbolics/expansion.jl")
include("Symbolics/pythagorean.jl")
include("Symbolics/fourier.jl")
include("Symbolics/harmonic_projection.jl")
include("Symbolics/drop_powers.jl")
include("Symbolics/linear_solve.jl")
include("DifferentialEquation.jl")
include("Variables.jl")
include("HarmonicVariable.jl")
include("HarmonicEquation.jl")

end
