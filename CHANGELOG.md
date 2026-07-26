# Changelog

All notable changes to QuestBase.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.0] - 2026-03-10

### Breaking

- Upgraded dependency bounds to SymbolicUtils.jl v4 and Symbolics.jl v7. Downstream packages must be compatible with these versions.
- Dropped support for SymbolicUtils.jl v3 and Symbolics.jl v6.

### Changed

- Replaced all uses of SymbolicUtils internals (`@compactified`, direct field access like `.base`, `.exp`, `.num`, `.den`, `.f`, `.arguments`, `.coeff`, `.name`) with public API (`isadd`, `ismul`, `isdiv`, `ispow`, `isterm`, `issym`, `operation`, `arguments`, `nameof`, `unwrap_const`).
- Replaced `Symbolics._toexpr` with the public `Symbolics.tosymbol` in `var_name`.
- Replaced `Symbolics.isterm`/`Symbolics.issym`/`Symbolics.is_derivative` with imports from SymbolicUtils/Symbolics rather than qualified calls.
- Removed dependency on `SymbolicUtils.Unityper` (removed in SymbolicUtils v4).
- Removed use of `frac_maketerm` in `add_div` (deleted in SymbolicUtils v4).
- Adapted `substitute_all` to handle Symbolics v7's change where `substitute` no longer recurses into `Differential` arguments; added `include_derivatives` keyword.
- Adapted `_parameters` in `HarmonicEquation` to handle `Symbolics.get_variables` returning `Set{BasicSymbolic}` instead of `Vector{Num}`, and to filter out derivative terms now returned by `get_variables` in v7.
- Adapted `rearrange!` in `DifferentialEquation` for v7's `get_variables` treating derivatives as variables.
- Adapted `trig_to_exp` to work at the `BasicSymbolic` level to avoid `Complex{Num}` issues with v7's `exp(im * Num)`.
- Adapted `exp_to_trig` to handle `Const`-wrapped zero arguments (`exp(Const(0)) = 1`).
- Adapted `_has_negative_coefficient` to unwrap `Const`-wrapped numeric coefficients in SymbolicUtils v4.
- Adapted `power_of`, `simplify_exp_products_mul`, and `trig_to_exp` to use `unwrap_const` for `Const`-wrapped numeric values in expression arguments.

### Performance

Benchmarks show significant improvements from SymbolicUtils v4 / Symbolics v7:

- **Fourier pipeline**: 3--5x faster (`fourier_cos_term`, `fourier_sin_term`, `trig_reduce`, `trig_to_exp`, `exp_to_trig`)
- **Power operations**: 3--4x faster (`max_power`, `drop_powers` with multiple variables)
- **Exponential simplification**: ~2x faster (`expand_exp_power`, `simplify_exp_products`)
- **Harmonic checks**: ~5x faster (`is_harmonic`)
- **Construction**: ~1.5x faster (`DifferentialEquation`)
- **Substitution**: 1.3--3x faster (`substitute_all`)
- Minor regressions on sub-millisecond operations (`expand_fraction`, `_apply_termwise` on `Div`, `drop_powers` with single variable)

### Added

- Benchmark suite in `benchmark/` for tracking performance across versions.
