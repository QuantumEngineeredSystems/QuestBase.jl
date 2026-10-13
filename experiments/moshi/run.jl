module MoshiExperiment

using BenchmarkTools: @benchmark, median
using Moshi.Match: @match
using QuestBase: QuestBase
using Symbolics: @variables, unwrap
using SymbolicUtils:
    SymbolicUtils,
    BasicSymbolic,
    AddMulVariant,
    arguments,
    isadd,
    ismul,
    isdiv,
    add_worker,
    mul_worker,
    vartype
using Test: @test, @testset

const BS = SymbolicUtils.BasicSymbolicImpl

# Change branch selection only. In particular, both implementations use the same
# bulk constructors; a sum/prod fold would introduce an unrelated quadratic cost.
function apply_match(f, x::BasicSymbolic)
    @match x begin
        BS.AddMul(; variant) => if variant == AddMulVariant.ADD
            add_worker(vartype(typeof(x)), map(f, arguments(x)))
        else
            mul_worker(vartype(typeof(x)), map(f, arguments(x)))
        end
        BS.Div(; num, den) => apply_match(f, num) / apply_match(f, den)
        _ => f(x)
    end
end

function simplify_match(x::BasicSymbolic)
    if isadd(x) || ismul(x) || isdiv(x)
        return apply_match(simplify_match, x)
    end
    return QuestBase.simplify_complex(x)
end
simplify_match(x) = QuestBase.simplify_complex(x)

function inputs()
    @variables a b c t x(t)
    return [
        "symbol" => unwrap(a),
        "power" => unwrap(a^3),
        "sum" => unwrap(a^2 + b * c + a * b^2 + c^3),
        "product" => unwrap(a * b * c * (a + b) * (b + c)),
        "fraction" => unwrap((a^2 + b * c + c^3) / (a * b + c^2)),
        "nested fraction" => unwrap((a + b / c) / (a * b + c)),
        "function" => unwrap(sin(a) + exp(b) + x),
        "wide sum" => unwrap(sum(a^i * b^(64-i) for i in 1:64)),
    ]
end

function run(; samples=100)
    cases = inputs()
    @testset "matching preserves production semantics" begin
        for (_, expression) in cases
            @test isequal(
                apply_match(identity, expression),
                QuestBase._apply_termwise(identity, expression),
            )
            @test isequal(
                apply_match(x -> x^2, expression),
                QuestBase._apply_termwise(x -> x^2, expression),
            )
            @test isequal(
                simplify_match(expression), QuestBase.simplify_complex(expression)
            )
        end
    end
    println("Julia ", VERSION, "; Julia threads=", Threads.nthreads())
    println(
        "QuestBase=", pkgversion(QuestBase), "; SymbolicUtils=", pkgversion(SymbolicUtils)
    )
    println("Steady-state medians; evals=1, samples=", samples, ", seconds≤1 per case")
    println(
        "| Workload | Predicates ns | Match ns | Predicates bytes | Match bytes | Predicates allocs | Match allocs |",
    )
    println("|:--|--:|--:|--:|--:|--:|--:|")
    rows = map(cases) do (name, expression)
        production = QuestBase.simplify_complex
        prototype = simplify_match
        # Warm both paths and check equality before collecting timings.
        production(expression)
        prototype(expression)
        before = @benchmark $production($expression) samples=samples evals=1 seconds=1
        after = @benchmark $prototype($expression) samples=samples evals=1 seconds=1
        p, m = median(before), median(after)
        println(
            "| ",
            name,
            " | ",
            round(p.time; digits=1),
            " | ",
            round(m.time; digits=1),
            " | ",
            p.memory,
            " | ",
            m.memory,
            " | ",
            p.allocs,
            " | ",
            m.allocs,
            " |",
        )
        return (; name, predicates=p, matching=m)
    end
    return rows
end

end

MoshiExperiment.run()
