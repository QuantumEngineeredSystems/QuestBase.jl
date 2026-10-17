# Symbolic Rewriting: Design Decisions

This document records what was tried when evaluating SymbolicUtils.jl's `@rule`-based
rewriting system as a replacement for QuestBase's manual tree traversal.

## Outcome

Rule-based rewriting was measured against manual traversal on QuestBase's typical
expressions (5 to 50 nodes). The rules lost in every case:

| Rewrite | `@rule` form | Slowdown |
|---|---|---|
| `exp(x)^n => exp(n*x)` | single rule + `Postwalk` | 1.9x |
| `exp(a)*exp(b) => exp(a+b)` | `@acrule` + `Fixpoint` + `Postwalk` | 6.8x |
| `(a+b)/c => a/c + b/c` | segment rule + `Postwalk` | 1.5x |
| exponential to trig, no-op input | custom node function + `Postwalk` | 2500x |

The exponential rewrites themselves have since been removed. Trigonometric reduction and
Fourier projection are now a single spectral pass (`src/Symbolics/harmonic_projection.jl`),
and linear solves run in a FLINT polynomial ring (`src/Symbolics/linear_solve.jl`).
Neither is a tree rewrite.

## Why It's Slower

1. **`Postwalk` visits every node**, even when only a few match. Manual code skips
   irrelevant subtrees.
2. **Pattern matching overhead.** `@rule` builds closures in continuation-passing style
   with `ImmutableDict` bindings. A check such as `ispow && isexp` is two boolean tests.
3. **`Fixpoint` re-traverses** the whole tree whenever rules interact.
4. **`@acrule` tries permutations** of the arguments up to `COMM_CHECKS_LIMIT`.

The rule system is designed for Symbolics.jl's simplifier, which applies hundreds of
rules to large trees. Its overhead is amortized there but dominates for small, targeted
transformations.

## Where `Postwalk` Fits

`collapse_pythagorean` uses `Postwalk(collapse_pythagorean_sum)`. The rewrite must see
every sum, and the node function is a hash lookup on squared trig factors, so the full
visit costs nothing extra.

## SymbolicUtils API Quick Reference

For QuestBase developers working with symbolic expressions.

### Expression Type Predicates
- `isadd(x)`, `ismul(x)`, `isdiv(x)`, `ispow(x)` — check structural type
- `isterm(x)` — true for Pow and Term nodes (compound expressions)
- `issym(x)` — true for leaf symbols
- `iscall(x)` — true for ALL compound nodes (more general than `isterm`)
- `is_literal_number(x)` — true for numeric literals
- `is_derivative(x)` — true for derivative expressions

### Expression Accessors
- `arguments(x)` — child expressions (order depends on auto-simplification)
- `sorted_arguments(x)` — children in deterministic order
- `operation(x)` — head operation (`+`, `*`, `sin`, `exp`, etc.)
- `unwrap_const(x)` — extract value from `Const`-wrapped numbers (SymbolicUtils v4)
- `unwrap(x)` / `wrap(x)` — convert between `Num` and `BasicSymbolic`

### Rule System
- `@rule LHS => RHS` — basic pattern rule with `~x` (slot) and `~~xs` (segment)
- `@acrule` — associative-commutative (tries permutations, expensive)
- `@ordered_acrule` — cheaper than `@acrule` (combinations instead of permutations)
- `@capture expr pattern` — one-off pattern match, injects bindings into scope
- `~x::predicate` — slot with predicate filter (e.g., `~x::is_literal_number`)
- `~!x` — default slot (defaults to 0 in `+`, 1 in `*`, 1 in `^` exponent)

### Rewriter Combinators
- `Chain((r1, r2, ...))` — apply in order (use **tuples** for `@generated` unrolling)
- `RestartedChain((r1, r2, ...))` — restart from beginning on any match
- `Fixpoint(rw)` — repeat until no change
- `Postwalk(rw)` / `Prewalk(rw)` — bottom-up / top-down tree traversal
- `PassThrough(rw)` — convert `nothing` return to identity (pass-through)
- `If(cond, rw)` / `IfElse(cond, rw1, rw2)` — conditional dispatch

### Performance Tips
- `COMM_CHECKS_LIMIT` (default 10): commutative matching disabled for 10+ args
- `Postwalk` uses copy-on-write: unchanged subtrees are not rebuilt
- Hash-consing makes `===` checks O(1) for `Fixpoint` termination
- Guard rule groups with `If(predicate, ...)` to skip irrelevant subtrees
