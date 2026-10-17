# Fourier projection of expressions polynomial in cos and sin of arguments linear in time.
#
# After `expand`, every addend of `x` reads `c * m * Π_j f_j(ν_j t + φ_j)^p_j / d`, where
# `m`, `d` and the frequencies `ν_j` and phases `φ_j` are free of `t`, `f_j` is `cos` or
# `sin`, and `p_j` is a positive integer. With `cos(a) = (e^{ia} + e^{-ia})/2` and
# `sin(a) = (-i e^{ia} + i e^{-ia})/2`, the trig part is a finite sum `Σ w e^{i(ν t + φ)}`
# with numeric weights `w`. The coefficients of `cos(ωt)` and `sin(ωt)` come from the terms
# with `ν = ±ω`, and each phase contributes `cos(φ)` and `sin(φ)` factors. Addends that
# share a trig part share its spectrum, which is therefore built once.
#
# `scan_projection` walks the expanded expression once, checks that it has this form and
# records its structure. Numeric constants stay symbolic during the scan, because
# SymbolicUtils stores their values untyped; the scan only determines the real type `R`
# that holds all of them exactly. `ProjectionProblem{R}` reads every constant into storage
# of type `R`, and `project` then runs on concrete types throughout.

"The addends of a sum, or `x` itself."
function summands_of(x::SymbolicNode)
    return @match x begin
        BSImpl.AddMul(; variant) && if variant === AddMulVariant.ADD
        end => arguments(x)
        _ => SymbolicArguments(SymbolicUtils.ArgsT{NUM_VARTYPE}((x,)))
    end
end

"The factors of a product, or `x` itself."
function factors_of(x::SymbolicNode)
    return @match x begin
        BSImpl.AddMul(; variant) && if variant === AddMulVariant.MUL
        end => arguments(x)
        _ => SymbolicArguments(SymbolicUtils.ArgsT{NUM_VARTYPE}((x,)))
    end
end

"The product of `factors`, which may be empty."
function monomial(factors::Vector{SymbolicNode})
    isempty(factors) && return SYMBOLIC_ONE
    length(factors) == 1 && return only(factors)
    return mul_worker(NUM_VARTYPE, copy(factors))
end

"The sum of `terms`, which may be empty."
function polynomial(terms::Vector{SymbolicNode})
    isempty(terms) && return SYMBOLIC_ZERO
    length(terms) == 1 && return only(terms)
    return add_worker(NUM_VARTYPE, terms)
end

@enum ProjectionKind ConstantPart CosinePart SinePart

"Symbolic monomials numbered in order of appearance."
struct AtomTable
    index::Dict{SymbolicNode,Int}
    atoms::Vector{SymbolicNode}
end
AtomTable() = AtomTable(Dict{SymbolicNode,Int}(), SymbolicNode[])

function atom_index!(table::AtomTable, atom::SymbolicNode)
    return get!(table.index, atom) do
        push!(table.atoms, atom)
        return length(table.atoms)
    end
end

"A linear combination of atoms, as `atom index => numeric constant` pairs."
const AtomCombination = Vector{Pair{Int,SymbolicNode}}

"The argument `ν t + φ` of a trig factor, split into frequency and phase."
struct TrigArgument
    frequency::AtomCombination
    phase::AtomCombination
end

"A trig factor `cos(ν t + φ)` or `sin(ν t + φ)`."
struct TrigFactor
    sine::Bool
    argument::TrigArgument
end

"A requested harmonic: the constant part, or the coefficient of `cos(ωt)` or `sin(ωt)`."
struct HarmonicTarget
    kind::ProjectionKind
    frequency::AtomCombination
end

"""
An addend `coefficient * Π rest * Π_k trig_k^p_k / denominator`, with its trig part stored
as `k => p_k` pairs sorted by `k`. The `t`-free factors `rest` are sorted by hash, so equal
products compare equal without being built as expressions.
"""
struct ScannedAddend
    coefficient::SymbolicNode
    rest::Vector{SymbolicNode}
    denominator::SymbolicNode
    signature::Vector{Pair{Int,Int}}
end

"""
The number type a projection is computed in, read off the constants of the expression:
`Rational{Int}` unless a constant needs a big integer (`big`) or is inexact (`inexact`).
"""
struct NumericKind
    big::Bool
    inexact::Bool
end

NumericKind(::Union{Int,Bool,Rational{Int}}) = NumericKind(false, false)
NumericKind(::Union{BigInt,Rational{BigInt}}) = NumericKind(true, false)
NumericKind(::Union{Float64,Float32,Irrational}) = NumericKind(false, true)
NumericKind(::BigFloat) = NumericKind(true, true)
function NumericKind(value::Complex)
    iszero(imag(value)) && return NumericKind(real(value))
    return throw(
        ArgumentError("Fourier projection requires real coefficients, found $value.")
    )
end

"The kind that holds the constants of both `a` and `b`."
Base.:|(a::NumericKind, b::NumericKind) = NumericKind(a.big | b.big, a.inexact | b.inexact)

"The expanded expression in the form `ProjectionProblem` consumes."
struct ScannedExpression
    addends::Vector{ScannedAddend}
    trig_factors::Vector{TrigFactor}
    targets::Vector{HarmonicTarget}
    frequencies::AtomTable
    phases::AtomTable
    numeric_kind::NumericKind
end

"State of `scan_projection`."
struct ProjectionScan
    time::SymbolicNode
    frequencies::AtomTable
    phases::AtomTable
    trig_index::Dict{SymbolicNode,Int}
    trig_factors::Vector{TrigFactor}
    time_dependence::Dict{SymbolicNode,Bool}
    numeric_kind::Base.RefValue{NumericKind}
end

function ProjectionScan(time::SymbolicNode)
    return ProjectionScan(
        time,
        AtomTable(),
        AtomTable(),
        Dict{SymbolicNode,Int}(),
        TrigFactor[],
        Dict{SymbolicNode,Bool}(),
        Ref(NumericKind(false, false)),
    )
end

"Return `true` if `x` contains the time variable of `scan`."
function depends_on_time(scan::ProjectionScan, x::SymbolicNode)
    return get!(scan.time_dependence, x) do
        isequal(x, scan.time) && return true
        iscall(x) || return false
        return any(argument -> depends_on_time(scan, argument), arguments(x))
    end
end

"Record the numeric kind of the constant `node` and return its product with `coefficient`."
function numeric_constant!(
    scan::ProjectionScan, coefficient::SymbolicNode, node::SymbolicNode
)
    scan.numeric_kind[] |= with_const(NumericKind, node)
    return isequal(coefficient, SYMBOLIC_ONE) ? node : coefficient * node
end

function non_polynomial(factor::SymbolicNode)
    return throw(ArgumentError("Cannot project the non-polynomial trig factor $factor."))
end

"The exponent `power` of `factor`, which must be a positive integer."
function positive_power(power::Integer, factor::SymbolicNode)
    return power > 0 ? Int(power) : non_polynomial(factor)
end
positive_power(power::Number, factor::SymbolicNode) = non_polynomial(factor)

"The exponent of a trig factor `f(ν t + φ)^p`, which must be a positive integer."
function trig_power(factor::SymbolicNode)
    exponent = last(base_exponent(factor))
    isconst(exponent) || non_polynomial(factor)
    return with_const(power -> positive_power(power, factor), exponent)
end

"""
Whether the trig factor `factor` enters the spectrum, which needs a positive integer power. A
time-free factor that fails this, such as `1/cos(ψ)`, is a plain coefficient instead.
"""
function is_spectral_power(factor::SymbolicNode)
    exponent = last(base_exponent(factor))
    return isconst(exponent) && with_const(is_positive_integer, exponent)
end
is_positive_integer(value::Integer) = value > 0
is_positive_integer(value::Number) = false

function nonlinear_argument(argument::SymbolicNode)
    return throw(ArgumentError("The trig argument $argument is not linear in time."))
end

"Parse `ν t + φ` into frequency and phase combinations, requiring linearity in `t`."
function scan_argument(scan::ProjectionScan, argument::SymbolicNode)
    frequency = AtomCombination()
    phase = AtomCombination()
    factors = SymbolicNode[]
    for term in expanded_terms(argument)
        depends_on_time(scan, term.denominator) && nonlinear_argument(argument)
        coefficient = SYMBOLIC_ONE
        time_power = 0
        empty!(factors)
        for factor in term.factors
            if isconst(factor)
                coefficient = numeric_constant!(scan, coefficient, factor)
            elseif isequal(factor, scan.time)
                time_power += 1
            elseif depends_on_time(scan, factor)
                nonlinear_argument(argument)
            else
                push!(factors, factor)
            end
        end
        time_power <= 1 || nonlinear_argument(argument)
        atom = term_expression(ExpandedTerm(factors, term.denominator))
        if isone(time_power)
            push!(frequency, atom_index!(scan.frequencies, atom) => coefficient)
        else
            push!(phase, atom_index!(scan.phases, atom) => coefficient)
        end
    end
    return TrigArgument(frequency, phase)
end

function trig_factor_index!(scan::ProjectionScan, trig::SymbolicNode)
    return get!(scan.trig_index, trig) do
        argument = scan_argument(scan, first(arguments(trig)))
        push!(scan.trig_factors, TrigFactor(operation(trig) === sin, argument))
        return length(scan.trig_factors)
    end
end

function scan_target(scan::ProjectionScan, ω, f)
    f === cos ||
        f === sin ||
        throw(ArgumentError("Harmonics are projected onto cos or sin, not $f."))
    isequal(ω, 0) &&
        return HarmonicTarget(f === cos ? ConstantPart : SinePart, AtomCombination())
    frequency = Num(ω)
    depends_on_time(scan, unwrap(frequency)) &&
        throw(ArgumentError("The harmonic frequency $ω must not depend on time."))
    argument = scan_argument(scan, unwrap(frequency * Num(scan.time)))
    return HarmonicTarget(f === cos ? CosinePart : SinePart, argument.frequency)
end

"Split one expanded term into its coefficient, `t`-free factors and trig part."
function scan_addend(scan::ProjectionScan, term::ExpandedTerm)
    depends_on_time(scan, term.denominator) && throw(
        ArgumentError("Cannot project: time enters the denominator $(term.denominator)."),
    )
    coefficient = SYMBOLIC_ONE
    rest = SymbolicNode[]
    signature = Pair{Int,Int}[]
    for factor in term.factors
        if isconst(factor)
            coefficient = numeric_constant!(scan, coefficient, factor)
            continue
        end
        base = first(base_exponent(factor))
        if is_trig_call(base) &&
            (is_spectral_power(factor) || depends_on_time(scan, factor))
            push!(signature, trig_factor_index!(scan, base) => trig_power(factor))
        elseif depends_on_time(scan, factor)
            throw(
                ArgumentError(
                    "Cannot project $factor: time enters other than through cos or sin."
                ),
            )
        else
            push!(rest, factor)
        end
    end
    sort!(rest; by=hash)
    return ScannedAddend(
        coefficient, rest, term.denominator, merge_powers!(sort!(signature))
    )
end

"The sorted signature `k => p_k` with the powers of a repeated trig factor `k` summed."
function merge_powers!(signature::Vector{Pair{Int,Int}})
    merged = 0
    for (index, power) in signature
        if merged > 0 && first(signature[merged]) == index
            signature[merged] = index => last(signature[merged]) + power
        else
            signature[merged += 1] = index => power
        end
    end
    return resize!(signature, merged)
end

"""
Record what `ProjectionProblem` needs from the expanded `terms`, or throw if they cannot be
projected.
"""
function scan_projection(terms::Vector{ExpandedTerm}, time, specs)
    scan = ProjectionScan(unwrap(Num(time)))
    targets = [scan_target(scan, ω, f) for (ω, f) in specs]
    addends = ScannedAddend[scan_addend(scan, term) for term in terms]
    return ScannedExpression(
        addends,
        scan.trig_factors,
        targets,
        scan.frequencies,
        scan.phases,
        scan.numeric_kind[],
    )
end

# Typed kernel. A spectrum maps `[frequency; phase]` coefficient vectors to weights.

const Spectrum{R} = Dict{Vector{R},Complex{R}}

"A target with its frequency as a dense coefficient vector."
struct DenseTarget{R}
    kind::ProjectionKind
    frequency::Vector{R}
end

"The scanned expression with every numeric constant read into type `R`."
struct ProjectionProblem{R}
    addends::Vector{ScannedAddend}
    coefficients::Vector{R}
    factor_spectra::Vector{Spectrum{R}}
    targets::Vector{DenseTarget{R}}
    phase_atoms::Vector{SymbolicNode}
    key_length::Int
end

function dense_combination(::Type{R}, combination::AtomCombination, n::Int) where {R}
    vector = zeros(R, n)
    for (index, coefficient) in combination
        vector[index] = with_const(R, coefficient)
    end
    return vector
end

function trig_spectrum(::Type{R}, factor::TrigFactor, nf::Int, np::Int) where {R}
    key = vcat(
        dense_combination(R, factor.argument.frequency, nf),
        dense_combination(R, factor.argument.phase, np),
    )
    half = one(R) / 2
    w₊, w₋ = if factor.sine
        Complex(zero(half), -half), Complex(zero(half), half)
    else
        Complex(half, zero(half)), Complex(half, zero(half))
    end
    return Spectrum{R}(key => w₊, -key => w₋)
end

function ProjectionProblem{R}(scanned::ScannedExpression) where {R}
    nf = length(scanned.frequencies.atoms)
    np = length(scanned.phases.atoms)
    return ProjectionProblem{R}(
        scanned.addends,
        R[with_const(R, addend.coefficient) for addend in scanned.addends],
        [trig_spectrum(R, factor, nf, np) for factor in scanned.trig_factors],
        [
            DenseTarget{R}(target.kind, dense_combination(R, target.frequency, nf)) for
            target in scanned.targets
        ],
        scanned.phases.atoms,
        nf + np,
    )
end

"""
    with_problem(f, scanned)

`f` applied to the `ProjectionProblem` of `scanned` in the number type its constants need.
Each branch names its type, so `f` is compiled for every one and the call stays inferable.
"""
function with_problem(f::F, scanned::ScannedExpression) where {F}
    kind = scanned.numeric_kind
    if kind.inexact
        kind.big && return f(ProjectionProblem{BigFloat}(scanned))
        return f(ProjectionProblem{Float64}(scanned))
    end
    kind.big && return f(ProjectionProblem{Rational{BigInt}}(scanned))
    return f(ProjectionProblem{Rational{Int}}(scanned))
end

"The constant `x`, with integer-valued rationals stored as integers."
symbolic_number(x::Real) = Const{NUM_VARTYPE}(x)
function symbolic_number(x::Rational)
    return Const{NUM_VARTYPE}(isone(denominator(x)) ? numerator(x) : x)
end

"`cosine * cos(φ) + sine * sin(φ)` for the phase `φ` numbered `phase`; phase 1 is zero."
struct PhaseContribution{R}
    phase::Int
    cosine::R
    sine::R
end

"""
Canonical phases, with `cos(φ)` and `sin(φ)` built once per phase. The atom
`SYMBOLIC_ONE`, if present, carries the numeric offset of a phase and is left out of the
canonical `order`.
"""
struct PhaseTable{R}
    index::Dict{Vector{R},Int}
    cosines::Vector{SymbolicNode}
    sines::Vector{SymbolicNode}
    atoms::Vector{SymbolicNode}
    order::Vector{Int}
end

function PhaseTable{R}(atoms::Vector{SymbolicNode}) where {R}
    symbolic = findall(atom -> !isequal(atom, SYMBOLIC_ONE), atoms)
    order = symbolic[sortperm(atoms[symbolic]; lt=SymbolicUtils.:<ₑ)]
    index = Dict(zeros(R, length(atoms)) => 1)
    return PhaseTable{R}(index, [SYMBOLIC_ONE], [SYMBOLIC_ZERO], atoms, order)
end

"The first nonzero symbolic coefficient of `phase` in the canonical order of its atoms."
function leading_coefficient(table::PhaseTable{R}, phase::Vector{R}) where {R}
    for i in table.order
        iszero(phase[i]) || return phase[i]
    end
    return zero(R)
end

"Build `cos(φ)` and `sin(φ)` for a canonical phase and return its number."
function add_phase!(table::PhaseTable{R}, phase::Vector{R}) where {R}
    # A purely numeric phase, whose only nonzero entry is the offset, is evaluated; a
    # symbolic one keeps its numeric offset exact.
    if iszero(leading_coefficient(table, phase))
        push!(table.cosines, Const{NUM_VARTYPE}(cos(sum(phase))))
        push!(table.sines, Const{NUM_VARTYPE}(sin(sum(phase))))
    else
        terms = [
            symbolic_number(c) * atom for (c, atom) in zip(phase, table.atoms) if !iszero(c)
        ]
        angle = polynomial(terms)
        push!(table.cosines, cos(angle))
        push!(table.sines, sin(angle))
    end
    return table.index[phase] = length(table.cosines)
end

"""
Number `phase` after fixing its sign so that its first nonzero symbolic coefficient is
positive. Return the number and whether the sign was flipped, which negates `sin(φ)` and
leaves `cos(φ)`.
"""
function phase_index!(table::PhaseTable{R}, phase::Vector{R}) where {R}
    flipped = leading_coefficient(table, phase) < 0
    canonical = flipped ? -phase : phase
    index = get(table.index, canonical, 0)
    iszero(index) && (index = add_phase!(table, canonical))
    return index, flipped
end

function convolve(a::Spectrum{R}, b::Spectrum{R}) where {R}
    product = Spectrum{R}()
    for (ka, wa) in a, (kb, wb) in b
        key = ka + kb
        product[key] = get(product, key, zero(Complex{R})) + wa * wb
    end
    return filter!(!iszero ∘ last, product)
end

"Whether the frequency part of `key` equals `+ω` and `-ω`."
function frequency_match(key::Vector{R}, ω::Vector{R}) where {R}
    plus = minus = true
    for i in eachindex(ω)
        plus &= key[i] == ω[i]
        minus &= key[i] == -ω[i]
    end
    return plus, minus
end

"""
Coefficients of one target in the spectrum `Σ w e^{i(ν t + φ)}`. With `w = a + ib`, a term
contributes `Re(w e^{iφ}) = a cos φ - b sin φ` to the constant part and to `cos(ωt)` when
`ν = ±ω`, and `∓Im(w e^{iφ}) = ∓(a sin φ + b cos φ)` to `sin(ωt)` when `ν = ±ω`.
"""
function project_spectrum(
    spectrum::Spectrum{R}, target::DenseTarget{R}, phases::PhaseTable{R}
) where {R}
    nf = length(target.frequency)
    totals = Dict{Int,Tuple{R,R}}()
    for (key, w) in spectrum
        plus, minus = frequency_match(key, target.frequency)
        plus | minus || continue
        a, b = real(w), imag(w)
        α, β = if target.kind != SinePart
            (a, -b)
        elseif plus && minus # sin(0 t) vanishes
            (zero(R), zero(R))
        else
            plus ? (-b, -a) : (b, a)
        end
        index, flipped = phase_index!(phases, key[(nf + 1):end])
        α₀, β₀ = get(totals, index, (zero(R), zero(R)))
        totals[index] = (α₀ + α, β₀ + (flipped ? -β : β))
    end
    contributions = PhaseContribution{R}[]
    for (index, (α, β)) in totals
        iszero(α) && iszero(β) || push!(contributions, PhaseContribution{R}(index, α, β))
    end
    return contributions
end

"The spectrum of the trig product `signature`."
function product_spectrum(
    problem::ProjectionProblem{R}, signature::Vector{Pair{Int,Int}}
) where {R}
    spectrum = Spectrum{R}(zeros(R, problem.key_length) => one(Complex{R}))
    for (factor, power) in signature, _ in 1:power
        spectrum = convolve(spectrum, problem.factor_spectra[factor])
    end
    return spectrum
end

"Projections of the trig product `signature` onto every target."
function signature_contributions(
    problem::ProjectionProblem{R}, signature::Vector{Pair{Int,Int}}, phases::PhaseTable{R}
) where {R}
    spectrum = product_spectrum(problem, signature)
    return [project_spectrum(spectrum, target, phases) for target in problem.targets]
end

"""
The numeric weight of every term `rest * trig` of one numerator. Like terms are summed here,
in exact arithmetic, so each surviving term is built as an expression only once.
"""
const TermWeights{R} = OrderedDict{Tuple{Vector{SymbolicNode},SymbolicNode},R}

"""
The terms of one projected coefficient, grouped by the denominator of the addend they came
from. Each group is summed before it is divided, so equal denominators are divided once and
their numerators merge, instead of leaving one fraction per term.
"""
const DenominatorGroups{R} = OrderedDict{SymbolicNode,TermWeights{R}}

function push_term!(
    groups::DenominatorGroups{R}, weight::R, addend::ScannedAddend, trig::SymbolicNode
) where {R}
    iszero(weight) && return groups
    weights = get!(TermWeights{R}, groups, addend.denominator)
    key = (addend.rest, trig)
    weights[key] = get(weights, key, zero(R)) + weight
    return groups
end

"The sum of every group's numerator over its denominator."
function grouped_sum(groups::DenominatorGroups{R}) where {R}
    fractions = SymbolicNode[]
    for (denominator, weights) in groups
        terms = SymbolicNode[
            mul_worker(NUM_VARTYPE, SymbolicNode[symbolic_number(weight); rest; trig]) for
            ((rest, trig), weight) in weights if !iszero(weight)
        ]
        isempty(terms) && continue
        numerator = polynomial(terms)
        push!(
            fractions,
            isequal(denominator, SYMBOLIC_ONE) ? numerator : numerator / denominator,
        )
    end
    return polynomial(fractions)
end

"The coefficient of every target of `problem`."
function project(problem::ProjectionProblem{R}) where {R}
    terms = [DenominatorGroups{R}() for _ in problem.targets]
    phases = PhaseTable{R}(problem.phase_atoms)
    cache = Dict{Vector{Pair{Int,Int}},Vector{Vector{PhaseContribution{R}}}}()
    for (addend, c) in zip(problem.addends, problem.coefficients)
        projected = get!(cache, addend.signature) do
            return signature_contributions(problem, addend.signature, phases)
        end
        for (target_terms, contributions) in zip(terms, projected), q in contributions
            push_term!(target_terms, c * q.cosine, addend, phases.cosines[q.phase])
            q.phase == 1 && continue
            push_term!(target_terms, c * q.sine, addend, phases.sines[q.phase])
        end
    end
    return SymbolicNode[grouped_sum(groups) for groups in terms]
end

"""
    fourier_terms(x, t, specs)

Coefficients of `f(ω*t)` in `x` for every `(ω, f)` pair in `specs`, where each `f` is
`cos` or `sin`. The pair `(0, cos)` selects the time-independent part.

`x` must be polynomial in `cos` and `sin` of arguments linear in `t`, with real
coefficients and `t` entering nowhere else; denominators must be free of `t`. Phases such
as `cos(ωt + ψ)` are allowed and produce `cos(ψ)` and `sin(ψ)` factors. Throws an
`ArgumentError` otherwise, since such inputs have no finite Fourier series.

The coefficients are read off the exact exponential spectrum of each expanded addend,
which is built once for every distinct trig product in `x`.
"""
function fourier_terms(x, t, specs)
    scanned = scan_projection(expanded_terms(unwrap(Num(x))), t, specs)
    return Num[Num(term) for term in with_problem(project, scanned)]
end

function fourier_terms(x::Complex{Num}, t, specs)
    isequal(x.im, 0) ||
        throw(ArgumentError("Fourier projection requires a real expression, found $x."))
    return fourier_terms(x.re, t, specs)
end

function fourier_terms(x::Equation, t, specs)
    lhs = fourier_terms(x.lhs, t, specs)
    rhs = fourier_terms(x.rhs, t, specs)
    return [Equation(l, r) for (l, r) in zip(lhs, rhs)]
end

"One harmonic of a Fourier series: `cosine * cos(frequency * t) + sine * sin(frequency * t)`."
struct FourierComponent
    frequency::Num
    cosine::Num
    sine::Num
end

"""
The first nonzero entry of the frequency vector `ν`, whose sign fixes the canonical sign of
the harmonic: `cos` is even and `sin` odd, so `-ν` is folded onto `ν`.
"""
function leading_entry(ν::AbstractVector{R}) where {R}
    for value in ν
        iszero(value) || return value
    end
    return zero(R)
end

"""
Accumulate every harmonic of `problem` into `cosines` and `sines`, keyed by the canonical
frequency vector.

A spectral term `w e^{i(ν t + φ)}` with `w = a + ib` contributes its real part
`(a cos φ - b sin φ) cos(νt) - (a sin φ + b cos φ) sin(νt)`. The spectrum of a real
expression pairs every term with its conjugate, so summing real parts over all terms gives
the series. Folding `ν < 0` onto `-ν` flips the sign of the sine part.
"""
function series!(
    cosines::Dict{Vector{R},DenominatorGroups{R}},
    sines::Dict{Vector{R},DenominatorGroups{R}},
    problem::ProjectionProblem{R},
    nf::Int,
) where {R}
    phases = PhaseTable{R}(problem.phase_atoms)
    cache = Dict{Vector{Pair{Int,Int}},Spectrum{R}}()
    for (addend, c) in zip(problem.addends, problem.coefficients)
        spectrum = get!(cache, addend.signature) do
            return product_spectrum(problem, addend.signature)
        end
        for (key, w) in spectrum
            ν = key[1:nf]
            folded = leading_entry(ν) < 0
            folded && (ν = -ν)
            index, flipped = phase_index!(phases, key[(nf + 1):end])
            a, b = real(w), imag(w)
            # sin(φ) of the stored canonical phase, which `flipped` negates
            s = flipped ? -one(R) : one(R)
            cosine = get!(DenominatorGroups{R}, cosines, ν)
            push_term!(cosine, c * a, addend, phases.cosines[index])
            index == 1 || push_term!(cosine, -c * b * s, addend, phases.sines[index])
            iszero(ν) && continue
            sign = folded ? one(R) : -one(R)
            sine = get!(DenominatorGroups{R}, sines, ν)
            push_term!(sine, sign * c * b, addend, phases.cosines[index])
            index == 1 || push_term!(sine, sign * c * a * s, addend, phases.sines[index])
        end
    end
    return cosines, sines
end

"The frequency `Σ ν_i atom_i` of a harmonic."
function frequency_expression(ν::AbstractVector{<:Real}, atoms::Vector{SymbolicNode})
    return polynomial(
        SymbolicNode[
            symbolic_number(value) * atom for
            (value, atom) in zip(ν, atoms) if !iszero(value)
        ],
    )
end

"""
    fourier_series(x, t)

The Fourier series of `x` in `t`, as one [`FourierComponent`](@ref) per harmonic present.
The time-independent part comes first, at frequency zero with a zero sine, followed by the
harmonics in increasing order of their frequency vectors. Harmonics whose coefficients
cancel are dropped.

`x` must satisfy the same conditions as for [`fourier_terms`](@ref).
"""
fourier_series(x, t) = expanded_fourier_series(expanded_terms(unwrap(Num(x))), t)

"[`fourier_series`](@ref) of the expanded `terms`."
function expanded_fourier_series(terms::Vector{ExpandedTerm}, t)
    scanned = scan_projection(terms, t, Tuple{Num,typeof(cos)}[])
    atoms = scanned.frequencies.atoms
    return with_problem(problem -> fourier_series(problem, atoms), scanned)
end

function fourier_series(
    problem::ProjectionProblem{R}, frequency_atoms::Vector{SymbolicNode}
) where {R}
    nf = length(frequency_atoms)
    cosines = Dict{Vector{R},DenominatorGroups{R}}()
    sines = Dict{Vector{R},DenominatorGroups{R}}()
    series!(cosines, sines, problem, nf)
    keys_sorted = sort!(collect(union(keys(cosines), keys(sines))))
    components = FourierComponent[]
    ordered = vcat(filter(iszero, keys_sorted), filter(!iszero, keys_sorted))
    for ν in ordered
        cosine = Num(grouped_sum(get(cosines, ν, DenominatorGroups{R}())))
        sine = Num(grouped_sum(get(sines, ν, DenominatorGroups{R}())))
        isequal(cosine, 0) && isequal(sine, 0) && continue
        push!(
            components,
            FourierComponent(Num(frequency_expression(ν, frequency_atoms)), cosine, sine),
        )
    end
    return components
end

"""
    fourier_cos_term(x, ω, t)

Coefficient of `cos(ω*t)` in `x`; `ω = 0` gives the time-independent part.
See [`fourier_terms`](@ref) for the admissible inputs.
"""
fourier_cos_term(x, ω, t) = only(fourier_terms(x, t, [(ω, cos)]))

"""
    fourier_sin_term(x, ω, t)

Coefficient of `sin(ω*t)` in `x`. See [`fourier_terms`](@ref) for the admissible inputs.
"""
fourier_sin_term(x, ω, t) = only(fourier_terms(x, t, [(ω, sin)]))
