# Symbolic Rewriting: Design Decisions

This document records what was tried when evaluating SymbolicUtils.jl's `@rule`-based
rewriting system as a replacement for QuestBase's manual tree traversal.

## What We Tried

We attempted to replace the manual `isadd`/`ismul`/`isdiv` dispatch in four functions
with `@rule` / `@acrule` + `Postwalk`:

### `expand_exp_power`

```julia
# Before: manual recursion (14 μs)
function expand_exp_power(expr::BasicSymbolic)
    if isadd(expr)
        sum(expand_exp_power(arg) for arg in arguments(expr))
    elseif ismul(expr)
        prod(expand_exp_power(arg) for arg in arguments(expr))
    elseif ispow(expr) && isexp(arguments(expr)[1])
        exp(arguments(arguments(expr)[1])[1] * arguments(expr)[2])
    else
        expr
    end
end

# After: single rule + Postwalk (28 μs, 1.9x slower)
const _expand_exp_power_rw = Postwalk(PassThrough(@rule exp(~x)^~n => exp(~x * ~n)))
expand_exp_power(expr::BasicSymbolic) = _expand_exp_power_rw(expr)
```

### `simplify_exp_products`

```julia
# Before: manual findall + index manipulation (10 μs)
# (see src/Symbolics/exponentials.jl)

# After: @acrule + Fixpoint + Postwalk (71 μs, 6.8x slower)
const _simplify_exp_products_rw = Fixpoint(
    Postwalk(
        PassThrough(
            Chain((
                @acrule(exp(~a) * exp(~b) => exp(~a + ~b)),
                @rule(exp(~x)^~n => exp(~x * ~n)),
                @rule(exp(~x::is_literal_number) => begin
                    iszero(unwrap_const(~x)) ? 1 : exp(~x)
                end),
            )),
        ),
    ),
)
```

The `Fixpoint` is needed because combining two exps can create a new `exp(x)^2` (via
hash-consing auto-simplification), which then needs the power rule, which may enable
further combining.

### `expand_fraction`

```julia
# Before: manual _apply_termwise dispatch (25 μs)

# After: single rule + Postwalk (37 μs, 1.5x slower)
const _expand_fraction_rw = Postwalk(
    PassThrough(@rule (+(~~xs)) / ~c => sum(x / ~c for x in ~~xs))
)
```

### `exp_to_trig`

```julia
# Before: manual _apply_termwise + inline conversion (0.01 μs on no-op input)

# After: Postwalk with custom node function (27 μs on no-op input, 2500x slower)
exp_to_trig(x::BasicSymbolic) = Postwalk(PassThrough(_exp_to_trig_node))(x)
```

The extreme regression on no-op input (expression with no `exp` nodes) happens because
`Postwalk` visits every node in the tree. The manual code sees `isadd` at the top level,
recurses into children, sees they're all `cos`/`sin` terms (not `exp`), and returns
immediately.

## Why It's Slower

For QuestBase's typical expressions (5-50 nodes):

1. **`Postwalk` visits every node** — even when only a few nodes match. Manual code
   skips irrelevant subtrees via `_apply_termwise`.
2. **Pattern matching overhead** — `@rule` builds closures with continuation-passing
   style and `ImmutableDict` bindings. For a simple `ispow && isexp` check, this is
   much heavier than two boolean tests.
3. **`Fixpoint` re-traverses** — When rules interact, `Fixpoint(Postwalk(...))` walks
   the full tree again. Manual code handles interactions in a single pass.
4. **`@acrule` tries permutations** — For `exp(a)*exp(b)`, it tries all argument
   orderings up to `COMM_CHECKS_LIMIT` (default 10).

The rule system is designed for Symbolics.jl's simplifier, which applies hundreds of
rules to large expression trees. The overhead is amortized there but dominates for
QuestBase's small, targeted transformations.

`Postwalk` is still used where it works well: `expand_all` and `add_div`, where the
rewriter needs to visit all nodes anyway.

## Function-by-Function Analysis

Every function in `src/Symbolics/` was analyzed for rule-based rewriting potential.
This section serves as a reference for future optimization attempts.

### Candidates That Were Tried (manual is faster for now)

| Function | `@rule` pattern | Slowdown | Notes |
|---|---|---|---|
| `expand_exp_power` | `exp(~x)^~n => exp(~x * ~n)` | 1.9x | Single pattern, `Postwalk` overhead dominates on small trees |
| `simplify_exp_products` | `@acrule exp(~a) * exp(~b) => exp(~a + ~b)` | 6.8x | Needs `Fixpoint` + 3 chained rules; n-ary `findall` is faster |
| `expand_fraction` | `(+(~~xs)) / ~c => sum(x/~c for x in ~~xs)` | 1.5x | Segment variable `~~xs` works but `_apply_termwise` is leaner |
| `exp_to_trig` | Custom node function via `PassThrough` | 2500x (no-op) | `_apply_termwise` skips non-exp subtrees; `Postwalk` visits all |

These could become faster than manual code if expression sizes grow significantly
(100+ nodes) or if SymbolicUtils optimizes the rewriter overhead in the future.

### Potentially Worth Revisiting

| Function | Possible `@rule` | Why it might help |
|---|---|---|
| `trig_to_exp` (both `Num` and `BasicSymbolic` versions) | `@rule sin(~x)^~n => im^n * ((exp(-im*~x) - exp(im*~x))/2)^n` and similar for `cos` | Currently builds a substitution dict manually by iterating `get_all_terms` and filtering trigs. A `Postwalk` with two rules would be more declarative. Worth trying if `trig_to_exp` becomes a bottleneck. |
| `simplify_complex` | `@rule` on `Const`-wrapped complex with zero imaginary | Currently uses `_apply_termwise` + `unwrap_const`. A single rule could handle the leaf case but traversal still needed. |
| `is_trig` | `@capture` macro | Could use `@capture expr (sin(~x) or cos(~x))` instead of manual `ispow` + `operation` check. Cleaner but unlikely faster. |

### Not Candidates (complex logic, not pattern-based)

| Function | File | Why rules don't fit |
|---|---|---|
| `get_independent` | `Symbolics_utils.jl` | Requires semantic `is_function()` check per argument — each term needs a variable-dependence query, not a structural pattern match |
| `exp_to_trig` (inner logic) | `fourier.jl` | Sign normalization uses alphabetic ordering of `sorted_arguments` and `_has_negative_coefficient` heuristics — too stateful for a rule |
| `drop_powers` | `drop_powers.jl` | Uses algebraic ε-substitution strategy: replaces `var → ε*var`, expands, removes high powers of ε. Not tree-pattern-based at all |
| `_fourier_term` | `fourier.jl` | Multi-step pipeline: `trig_reduce` → `get_independent` → complex simplification. Orchestration logic, not a rewrite |
| `_get_all_terms` | `Symbolics_utils.jl` | Collects terms into a flat list — accumulation, not transformation |
| `power_of` / `max_power` | `drop_powers.jl` | Query functions returning integers, not expression rewrites |
| `count_derivatives` | `Symbolics_utils.jl` | Single-node field access (`op.order`), not a traversal |
| `is_harmonic` | `Symbolics_utils.jl` | Predicate combining `get_all_terms`, `is_trig`, and `max_power` — orchestration |

### Already Using `Postwalk` Effectively

| Function | How | Why it works |
|---|---|---|
| `expand_all` | `Postwalk(expand_exp_power)` after `SymbolicUtils.expand` | Needs to visit all nodes; `expand_exp_power` is a simple leaf check |
| `add_div` | `Postwalk(add_with_div)` | Uses SymbolicUtils-provided `add_with_div` rewriter |

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
