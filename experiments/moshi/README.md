# Moshi branch-selection experiment

This experiment compares predicate dispatch with Moshi matching on SymbolicUtils 4 nodes. It does not change QuestBase runtime code or package dependencies. The prototype retains the production bulk constructors and recursion so that the comparison does not accidentally measure a quadratic `sum`/`prod` reconstruction instead of branch selection.

Run on Julia 1.13 from the repository root:

```sh
julia --project=experiments/moshi -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments/moshi experiments/moshi/run.jl
```

The environment uses the current checkout through `[sources]`. Keep the generated, ignored Manifest when reproducing a measurement; the default-branch sources elsewhere in the package can advance. Record `Pkg.status(; mode=Pkg.PKGMODE_MANIFEST)` alongside results when sharing them.

Before benchmarking, the runner checks exact symbolic equality with production for symbols, powers, sums, products, nested fractions, functions, and a 64-term polynomial. It then warms each implementation and reports median time, bytes, and allocation counts for 100 single-evaluation samples. Package loading and compilation are excluded from these steady-state timings.

The benchmark deliberately reaches private QuestBase and SymbolicUtils construction helpers. Their shape is an experimental dependency, not a supported package API. These measurements assess one traversal workload; they do not establish that all symbolic rewriting should use matching. A runtime change would need broader workloads and inference checks first.
