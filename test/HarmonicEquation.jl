# test/test_HarmonicEquation.jl
using Test
using Symbolics
using SymbolicUtils
using QuestBase
using QuestBase:
    DifferentialEquation,
    HarmonicVariable,
    HarmonicEquation,
    _parameters,
    get_variables,
    get_independent_variables,
    rearrange_standard,
    substitute_all,
    is_rearranged,
    _remove_brackets,
    declare_variables,
    dummy_symbolic_Jacobian,
    @eqtest,
    get_all_terms,
    get_independent

# Setup common test variables
@variables t, T
@variables x(t) y(t) u(T) v(T) w(T)
D = Differential(T)

# Create simple test equation
eq1 = D(u) ~ u + v
eq2 = D(v) ~ -u + v
nat_eq = DifferentialEquation([eq1, eq2], [x, y])

# Create test harmonic variables
hv1 = HarmonicVariable(u, "test", "u", Num(1.0), x)
hv2 = HarmonicVariable(v, "test", "v", Num(1.0), y)

# Test constructor
@testset "Construction" begin
    heq1 = HarmonicEquation([eq1, eq2], [hv1, hv2], nat_eq)
    heq2 = HarmonicEquation([eq1, eq2], [hv1, hv2], Num[], nat_eq)
    for heq in [heq1, heq2]
        heq = heq1
        @test heq.equations == [eq1, eq2]
        @test heq.variables == [hv1, hv2]
        @test heq.natural_equation == nat_eq
        @test heq.parameters == Num[]
        @test heq.jacobian isa Matrix{Num}
    end

    heq3 = HarmonicEquation([eq1, eq2], [hv1, hv2], Num[], Num[1 1; 1 1])
    @test isempty(heq3.natural_equation.harmonics)
end

@testset "Parameter handling" begin
    parvar = @variables p q
    eq_with_params = D(u) ~ p * v + q * v
    heq = HarmonicEquation([eq_with_params, eq2], [hv1, hv2], nat_eq)
    params = _parameters(heq)
    @eqtest params == parvar
end

@testset "Variable handling" begin
    heq = HarmonicEquation([eq1, eq2], [hv1, hv2], nat_eq)
    vars = get_variables(heq)
    @test length(vars) == 2
    @eqtest [T] == get_independent_variables(heq)
end

@testset "Equation manipulation" begin
    heq = HarmonicEquation([eq1, eq2], [hv1, hv2], nat_eq)

    # Test rearrange
    rearranged = rearrange_standard(heq)
    @test is_rearranged(rearranged)

    # Test substitute_all
    @variables a
    rules = Dict(u => a)
    subbed = substitute_all(heq, rules)
    @test !isequal(subbed.equations, heq.equations)
    @test !isequal(subbed.variables, heq.variables)
    @eqtest subbed.variables[1].symbol == a
end

@testset "Utility functions" begin
    heq = HarmonicEquation([eq1, eq2], [hv1, hv2], nat_eq)

    # Test _remove_brackets
    no_brackets = _remove_brackets(heq)

    list = get_all_terms.(Symbolics.expand_derivatives.(no_brackets))
    list = unique(filter(x -> !(x isa Real), Symbolics.unwrap.(reduce(vcat, list))))
    @test all(map(x -> !hasproperty(x, :arguments), list))

    # Test declare_variables
    declared = declare_variables(heq)
    list = get_all_terms.(Symbolics.expand_derivatives.(declared))
    list = unique(filter(x -> !(x isa Real), Symbolics.unwrap.(reduce(vcat, list))))
    @test all(map(x -> !hasproperty(x, :arguments), list))
end

@testset "is_rearranged MF path" begin
    # No derivative terms on either side => treated as arranged by construction (MF_bool).
    eq_alg1 = u ~ u + v
    eq_alg2 = v ~ u
    heq_alg = HarmonicEquation([eq_alg1, eq_alg2], [hv1, hv2], nat_eq)
    @test is_rearranged(heq_alg)
end

@testset "Hopf variable filtered from ansatz" begin
    hv_hopf = HarmonicVariable(w, "test", "Hopf", Num(1.0), x)
    heq_hopf = HarmonicEquation([eq1, eq2], [hv1, hv2, hv_hopf], nat_eq)
    ans = QuestBase._show_ansatz(heq_hopf)
    @test !occursin(string(hv_hopf.symbol), ans)
end
