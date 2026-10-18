using QuestBase
using Test
using JET
using Symbolics: Symbolics, Differential, Num, @variables, unwrap
using QuestBase:
    DifferentialEquation,
    HarmonicVariable,
    CosineQuadrature,
    SineQuadrature,
    create_harmonic_variable,
    kind_label,
    show_ansatz,
    HarmonicEquation,
    add_harmonic!,
    collapse_pythagorean,
    count_derivatives,
    d,
    declare_variable,
    declare_variables,
    drop_powers,
    dummy_symbolic_Jacobian,
    expand_fraction,
    fourier_cos_term,
    fourier_series,
    fourier_sin_term,
    fourier_terms,
    fraction_free_linear_solve,
    get_all_terms,
    get_independent,
    get_independent_variables,
    get_variables,
    is_harmonic,
    is_rearranged,
    max_power,
    nonpolynomial_terms,
    rearrange_standard,
    reduce_denominator,
    simplify_complex,
    solve_linear_system,
    substitute_all,
    trig_reduce,
    var_name

const TARGETS = (QuestBase,)

report_text(report) = sprint(show, MIME("text/plain"), report)

# Prints the offending reports on failure: a bare `@test isempty(...)` shows only the count.
function test_no_reports(report)
    reports = map(report_text, JET.get_reports(report))
    isempty(reports) || @info "unexpected JET reports" reports
    @test isempty(reports)
    return nothing
end

# Every entry must be free of errors (`report_call`) and of runtime dispatch (`report_opt`).
function test_entries(entries)
    for (name, thunk) in entries
        @testset "$name" begin
            test_no_reports(JET.@report_call target_modules = TARGETS thunk())
            test_no_reports(JET.@report_opt target_modules = TARGETS thunk())
        end
    end
    return nothing
end

@testset "JET report_package" begin
    test_no_reports(JET.report_package(QuestBase; target_modules=TARGETS))
end

@testset "JET per feature" begin
    @variables t T ω ω0 γ α F a b x y
    @variables X(t) U(T) V(T)
    D = Differential(T)

    duffing = d(X, t, 2) + γ * d(X, t) + ω0^2 * X + α * X^3 ~ F * cos(ω * t)
    diff_eq = DifferentialEquation(duffing, X)
    add_harmonic!(diff_eq, X, ω)
    u = HarmonicVariable(U, "test", CosineQuadrature, Num(1.0), X)
    v = HarmonicVariable(V, "test", SineQuadrature, Num(1.0), X)
    harmonic_eq = HarmonicEquation([D(U) + U ~ V, D(V) ~ -U + a * V], [u, v], diff_eq)

    trig_polynomial = (a * cos(ω * t) + b * sin(2ω * t))^3
    linear_system = [a * x + b * y ~ 1, b * x - a * y ~ 0]

    @testset "Construction" begin
        test_entries([
            "DifferentialEquation" => () -> DifferentialEquation(duffing, X),
            "add_harmonic!" => () -> add_harmonic!(deepcopy(diff_eq), X, ω),
            "get_independent_variables" => () -> get_independent_variables(diff_eq),
            "is_rearranged" => () -> is_rearranged(harmonic_eq),
            "get_variables" => () -> get_variables(harmonic_eq),
            "declare_variables" => () -> declare_variables(harmonic_eq),
            "declare_variable" => () -> declare_variable("z", t),
            "dummy_symbolic_Jacobian" => () -> dummy_symbolic_Jacobian(2),
            "create_harmonic_variable" =>
                () -> create_harmonic_variable(X, ω, t, CosineQuadrature; new_symbol="w"),
            "kind_label" => () -> kind_label(SineQuadrature),
            "show_ansatz" => () -> show_ansatz(harmonic_eq),
        ])
    end

    @testset "Symbolic utilities" begin
        test_entries([
            "substitute_all" => () -> substitute_all(harmonic_eq, Dict(a => 1)),
            "drop_powers" => () -> drop_powers(trig_polynomial, [a, b], 2),
            "max_power" => () -> max_power(trig_polynomial, a),
            "collapse_pythagorean" =>
                () -> collapse_pythagorean(a * cos(x)^2 + a * sin(x)^2),
            "simplify_complex" => () -> simplify_complex(unwrap(a + 0im)),
            "get_independent" => () -> get_independent(a + b * cos(ω * t), t),
            "get_all_terms" => () -> get_all_terms(trig_polynomial),
            "is_harmonic" => () -> is_harmonic(trig_polynomial, t),
            "count_derivatives" => () -> count_derivatives(d(X, t, 2)),
            "expand_fraction" => () -> expand_fraction((a + b) / x),
            "d" => () -> d(X^2, t),
            "var_name" => () -> var_name(U),
            "nonpolynomial_terms" => () -> nonpolynomial_terms(diff_eq),
        ])
    end

    @testset "Fourier projection" begin
        test_entries([
            "trig_reduce" => () -> trig_reduce(trig_polynomial),
            "reduce_denominator" => () -> reduce_denominator(a / (cos(x)^2 + sin(x)^2)),
            "fourier_series" => () -> fourier_series(trig_polynomial, t),
            "fourier_terms" =>
                () -> fourier_terms(trig_polynomial, t, [(ω, cos), (ω, sin)]),
            "fourier_cos_term" => () -> fourier_cos_term(trig_polynomial, ω, t),
            "fourier_sin_term" => () -> fourier_sin_term(trig_polynomial, ω, t),
        ])
    end

    @testset "Linear solves and rearrangement" begin
        test_entries([
            "fraction_free_linear_solve" =>
                () -> fraction_free_linear_solve(linear_system, [x, y]),
            "solve_linear_system" => () -> solve_linear_system(linear_system, [x, y]),
            "rearrange_standard DifferentialEquation" => () -> rearrange_standard(diff_eq),
            "rearrange_standard HarmonicEquation" => () -> rearrange_standard(harmonic_eq),
        ])
    end
end
