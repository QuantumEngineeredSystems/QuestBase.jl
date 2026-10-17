"""
Raised when Bareiss elimination is abandoned because its intermediate polynomials grew past
what is workable.

Elimination over the polynomial ring is exact, so a size limit is the only way the solve gives
up. A dense system in many symbolic atoms, such as the harmonic equations of a van der Pol
oscillator expanded in three harmonics, grows combinatorially even though nothing has gone
wrong. Symbolic LU is no refuge on exactly these systems, since it grows nested fractions
instead. This propagates to the caller, letting one that has a cheaper route
(`add_jacobian!` evaluates the Jacobian implicitly) take it.
"""
struct BareissTooLarge <: Exception
    message::String
end

function Base.showerror(io::IO, error::BareissTooLarge)
    return print(io, error.message)
end

"""
Largest number of monomials an intermediate Bareiss entry may reach before the solve gives
up with a [`BareissTooLarge`](@ref).

The systems this solve handles well stay far below the budget: 10 monomials for a Duffing
oscillator, 28 for a parametron, 344 for a van der Pol with two harmonics.
"""
const BAREISS_TERM_BUDGET = 2_000

"""
Ceiling on the monomial count a single polynomial product may be predicted to reach.

The size checks on finished entries are not enough on their own: one multiplication of two
wide polynomials can exhaust memory before its result can be measured. The number of terms
in a product is at most the product of the factors' term counts, so that bound is checked
before multiplying. It is deliberately far looser than [`BAREISS_TERM_BUDGET`](@ref),
because the bound is pessimistic and systems that solve fine reach it: a van der Pol with
two harmonics multiplies 344-term entries, predicting over a hundred thousand terms for a
result that is much smaller.
"""
const BAREISS_PRODUCT_BUDGET = 4_000_000

"""
    fraction_free_linear_solve(equations, variables)

Solve a square symbolic linear system using pivoted Bareiss–Jordan elimination.

The coefficients are read into a FLINT polynomial ring over ``ℚ`` (see [`RingAtoms`](@ref)).
Bareiss elimination then performs only exact polynomial divisions and produces one
determinant denominator per solution, instead of the nested fractions produced by symbolic
LU. Every division is exact, since the ring is an integral domain, so the solve succeeds on
any nonsingular system and only gives up on size (see [`BareissTooLarge`](@ref)).

The system is first split into independent subsystems (see [`system_blocks`](@ref)) and each
is eliminated on its own: for a harmonic ansatz that is `n` separate `2 × 2` solves instead
of one `2n × 2n`.
"""
function fraction_free_linear_solve(equations, variables)
    coefficients, offsets, islinear = Symbolics.linear_expansion(equations, variables)
    islinear || throw(ArgumentError("the equation system is not linear in the variables"))

    rows, columns = size(coefficients)
    rows == columns == length(offsets) ||
        throw(DimensionMismatch("fraction-free linear solve requires a square system"))

    right_hand_side = -offsets
    solution = Vector{Num}(undef, rows)
    for (block_rows, block_columns) in system_blocks(coefficients)
        solution[block_columns] = bareiss_solve(
            coefficients[block_rows, block_columns], right_hand_side[block_rows]
        )
    end
    return solution
end

"""
    solve_linear_system(equations, variables)

Solve the linear system `equations` for `variables`, the kernel behind every `rearrange!`.

The solve is [`fraction_free_linear_solve`](@ref). SymbolicUtils 4 leaves
`cos(x)^2 + sin(x)^2` standing, and these solutions carry the coefficient determinant in
their denominator, so the Pythagorean identity is collapsed in both numerator and
denominator. Without it the expressions grow every time a system is rearranged.
"""
function solve_linear_system(equations, variables)
    solution = fraction_free_linear_solve(equations, variables)
    return reduce_denominator.(collapse_pythagorean.(solution))
end

"Is `entry` zero on its face, so the row and column it joins are unrelated?"
is_structural_zero(entry::Num) = is_structural_zero(unwrap(entry))
is_structural_zero(entry::SymbolicNode) = isequal(entry, SYMBOLIC_ZERO)

"""
Split a system into the independent subsystems its sparsity pattern allows.

The coefficient matrix of an `n`-harmonic ansatz is `n` uncoupled `2×2` blocks. Solving it
whole is worse than wasteful: a Bareiss entry carries the determinant of everything
eliminated so far, so the solution comes out over the determinant of the full matrix rather
than the one block determinant that survives.

Blocks are the connected components of the bipartite graph joining row `r` to column `c`
where the coefficient is not [`is_structural_zero`](@ref). Returns `(rows, columns)` index
pairs. A system with one component, or with an unbalanced one (structurally singular, left
to the elimination to reject), is returned whole as a single block.
"""
function system_blocks(coefficients::AbstractMatrix)
    n = size(coefficients, 1)
    whole = [(collect(1:n), collect(1:n))]
    parent = collect(1:(2n))
    for column in 1:n, row in 1:n
        is_structural_zero(coefficients[row, column]) && continue
        row_root = union_find_root!(parent, row)
        column_root = union_find_root!(parent, n + column)
        row_root == column_root || (parent[row_root] = column_root)
    end

    block_of_root = Dict{Int,Int}()
    block_rows = Vector{Int}[]
    for row in 1:n
        block = get!(block_of_root, union_find_root!(parent, row)) do
            push!(block_rows, Int[])
            return length(block_rows)
        end
        push!(block_rows[block], row)
    end
    length(block_rows) > 1 || return whole

    block_columns = [Int[] for _ in block_rows]
    for column in 1:n
        block = get(block_of_root, union_find_root!(parent, n + column), 0)
        iszero(block) && return whole
        push!(block_columns[block], column)
    end
    all(length.(block_rows) .== length.(block_columns)) || return whole

    return collect(zip(block_rows, block_columns))
end

"Root of `index` in the union-find forest `parent`, compressing the path on the way."
function union_find_root!(parent::Vector{Int}, index::Int)
    root = index
    while parent[root] != root
        root = parent[root]
    end
    while parent[index] != root
        parent[index], index = root, parent[index]
    end
    return root
end

"""
The symbolic atoms of a linear system, the generators of the polynomial ring it is solved in.

Every subexpression that is not a sum, a product, a non-negative integer power or an exact
rational constant becomes an atom: symbols, `cos(ωt)`, `exp(x)`, `1/x`, but also floating
point, irrational and complex constants. Keeping inexact constants symbolic leaves the ring
over ``ℚ`` and every division in it exact; they are put back when the solution is read out.
"""
struct RingAtoms
    index::Dict{SymbolicNode,Int}
    atoms::Vector{SymbolicNode}
end
RingAtoms() = RingAtoms(Dict{SymbolicNode,Int}(), SymbolicNode[])

"Is the constant `value` an exact rational the coefficient field holds?"
is_ring_coefficient(value::Union{Integer,Rational}) = true
is_ring_coefficient(value::Number) = false

"Is the constant `value` an exponent the polynomial ring can raise to?"
is_ring_exponent(value::Integer) = value >= 0
is_ring_exponent(value::Number) = false

"Is the node `exponent` a constant the polynomial ring can raise to?"
function is_ring_exponent(exponent::SymbolicNode)
    return isconst(exponent) && with_const(is_ring_exponent, exponent)
end

"Is `x` a constant the ring represents, rather than an atom?"
is_ring_constant(x::SymbolicNode) = isconst(x) && with_const(is_ring_coefficient, x)

"Record the atoms of `x` in `atoms`, in order of first appearance."
function collect_atoms!(atoms::RingAtoms, x::SymbolicNode)
    is_ring_constant(x) && return atoms
    @match x begin
        BSImpl.AddMul(;) => for argument in arguments(x)
            collect_atoms!(atoms, argument)
        end
        BSImpl.Term(; f, args) && if f === (^) && is_ring_exponent(args[2])
        end => collect_atoms!(atoms, args[1])
        _ => get!(atoms.index, x) do
            push!(atoms.atoms, x)
            return length(atoms.atoms)
        end
    end
    return atoms
end
collect_atoms!(atoms::RingAtoms, x::Num) = collect_atoms!(atoms, unwrap(x))

"The exact constant `x` as an element of the coefficient field ``ℚ``."
ring_coefficient(x::SymbolicNode) = with_exact(QQ, x)

"The non-negative integer constant `x` as an exponent of the polynomial ring."
ring_exponent(x::SymbolicNode) = with_exact(Int, x)

"""
The polynomial ring a system is solved in, with one generator per atom of `atoms`.

The degree-lexicographic order makes the leading coefficient, which fixes the sign of each
solution's denominator, independent of how the atoms were numbered.
"""
function solution_ring(atoms::RingAtoms)
    ring, generators = Nemo.polynomial_ring(
        QQ, length(atoms.atoms); internal_ordering=:deglex, cached=false
    )
    return ring, generators
end

"`x` as an element of `ring`, whose generators are those of `atoms`."
function to_ring(
    ring::QQMPolyRing,
    generators::Vector{QQMPolyRingElem},
    atoms::RingAtoms,
    x::SymbolicNode,
)
    is_ring_constant(x) && return ring(ring_coefficient(x))
    return @match x begin
        BSImpl.AddMul(; variant) => if variant === AddMulVariant.ADD
            total = zero(ring)
            for argument in arguments(x)
                total += to_ring(ring, generators, atoms, argument)
            end
            total
        else
            product = one(ring)
            for argument in arguments(x)
                product *= to_ring(ring, generators, atoms, argument)
            end
            product
        end
        BSImpl.Term(; f, args) && if f === (^) && is_ring_exponent(args[2])
        end => to_ring(ring, generators, atoms, args[1])^ring_exponent(args[2])
        _ => generators[atoms.index[x]]
    end
end
function to_ring(
    ring::QQMPolyRing, generators::Vector{QQMPolyRingElem}, atoms::RingAtoms, x::Num
)
    return to_ring(ring, generators, atoms, unwrap(x))
end

"The rational `value` as a symbolic constant, in the narrowest machine type that holds it."
function symbolic_number(value::QQFieldElem)
    num, den = Nemo.numerator(value), Nemo.denominator(value)
    if isone(den)
        Nemo.fits(Int, num) && return Const{NUM_VARTYPE}(Int(num))
        return Const{NUM_VARTYPE}(BigInt(num))
    end
    if Nemo.fits(Int, num) && Nemo.fits(Int, den)
        return Const{NUM_VARTYPE}(Int(num) // Int(den))
    end
    return Const{NUM_VARTYPE}(BigInt(num) // BigInt(den))
end

"""
Rebuild a symbolic expression from a polynomial of the solution ring.

A solved entry carries thousands of monomials, so the sum and the per-monomial product go to
`add_worker`/`mul_worker` whole rather than through `+=`/`*=`. See [`_apply_termwise`](@ref).
"""
function from_ring(polynomial::QQMPolyRingElem, atoms::RingAtoms)
    addends = SymbolicNode[]
    factors = SymbolicNode[]
    for (coefficient, exponents) in
        zip(Nemo.coefficients(polynomial), Nemo.exponent_vectors(polynomial))
        empty!(factors)
        push!(factors, symbolic_number(coefficient))
        for (atom, exponent) in zip(atoms.atoms, exponents)
            iszero(exponent) || push!(factors, atom^exponent)
        end
        push!(
            addends,
            length(factors) == 1 ? only(factors) : mul_worker(NUM_VARTYPE, copy(factors)),
        )
    end
    isempty(addends) && return Num(SYMBOLIC_ZERO)
    return Num(length(addends) == 1 ? only(addends) : add_worker(NUM_VARTYPE, addends))
end

"Solve one block `coefficients * x = right_hand_side` by Bareiss–Jordan elimination."
function bareiss_solve(coefficients::AbstractMatrix, right_hand_side::AbstractVector)
    n = size(coefficients, 1)
    atoms = RingAtoms()
    foreach(entry -> collect_atoms!(atoms, entry), coefficients)
    foreach(entry -> collect_atoms!(atoms, entry), right_hand_side)
    ring, generators = solution_ring(atoms)

    augmented = Matrix{QQMPolyRingElem}(undef, n, n + 1)
    for column in 1:n, row in 1:n
        augmented[row, column] = to_ring(ring, generators, atoms, coefficients[row, column])
    end
    for row in 1:n
        augmented[row, n + 1] = to_ring(ring, generators, atoms, right_hand_side[row])
    end
    bareiss_jordan!(augmented)

    solution = Vector{Num}(undef, n)
    for row in 1:n
        numerator = augmented[row, n + 1]
        denominator = augmented[row, row]
        if Nemo.leading_coefficient(denominator) < 0
            numerator, denominator = -numerator, -denominator
        end
        solution[row] = from_ring(numerator, atoms) / from_ring(denominator, atoms)
    end
    return solution
end

"""
Bring the augmented matrix `[A | b]` to diagonal form in place by fraction-free
Gauss–Jordan elimination.

Afterwards every diagonal entry is `det(A)` and the last column holds `det(A) * x`. Each step
divides by the previous pivot, which Sylvester's identity guarantees to be exact, so the
entries stay polynomials. Row `pivot_index` is only read during its step, so the update runs
in place without a copy.
"""
function bareiss_jordan!(augmented::Matrix{QQMPolyRingElem})
    n = size(augmented, 1)
    previous_pivot = one(augmented[1, 1])
    for pivot_index in 1:n
        pivot_row = bareiss_pivot_row(augmented, pivot_index)
        iszero(augmented[pivot_row, pivot_index]) &&
            throw(LinearAlgebra.SingularException(pivot_index))
        if pivot_row != pivot_index
            for column in axes(augmented, 2)
                augmented[pivot_index, column], augmented[pivot_row, column] = augmented[
                    pivot_row, column
                ],
                augmented[pivot_index, column]
            end
        end

        pivot = augmented[pivot_index, pivot_index]
        for row in 1:n
            row == pivot_index && continue
            factor = augmented[row, pivot_index]
            for column in axes(augmented, 2)
                column == pivot_index && continue
                numerator =
                    bareiss_product(pivot, augmented[row, column], pivot_index, n) -
                    bareiss_product(factor, augmented[pivot_index, column], pivot_index, n)
                entry = Nemo.divexact(numerator, previous_pivot)
                length(entry) > BAREISS_TERM_BUDGET && throw(
                    BareissTooLarge(
                        "Bareiss entry reached $(length(entry)) terms at pivot " *
                        "$(pivot_index) of $(n), over the $(BAREISS_TERM_BUDGET) term budget",
                    ),
                )
                augmented[row, column] = entry
            end
            augmented[row, pivot_index] = zero(pivot)
        end
        previous_pivot = pivot
    end
    return augmented
end

"""
Row at or below `pivot_index` whose entry in column `pivot_index` has the fewest terms, so the
products that follow stay small. Every candidate is zero when the matrix is singular, and
the row returned then holds a zero.
"""
function bareiss_pivot_row(augmented::Matrix{QQMPolyRingElem}, pivot_index::Int)
    best_row = pivot_index
    best_size = typemax(Int)
    for row in pivot_index:size(augmented, 1)
        candidate = augmented[row, pivot_index]
        iszero(candidate) && continue
        if length(candidate) < best_size
            best_row = row
            best_size = length(candidate)
        end
    end
    return best_row
end

"Multiply two Bareiss entries, refusing products predicted to be unmanageably wide."
function bareiss_product(
    left::QQMPolyRingElem, right::QQMPolyRingElem, pivot_index::Int, n::Int
)
    predicted = length(left) * length(right)
    predicted > BAREISS_PRODUCT_BUDGET && throw(
        BareissTooLarge(
            "Bareiss product predicted at $(predicted) terms at pivot $(pivot_index) of " *
            "$(n), over the $(BAREISS_PRODUCT_BUDGET) term product budget",
        ),
    )
    return left * right
end
