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
