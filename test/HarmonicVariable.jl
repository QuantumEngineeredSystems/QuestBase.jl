# test/test_HarmonicVariable.jl
using Test
using Symbolics
using QuestBase:
    declare_variable,
    declare_variables,
    substitute_all,
    var_name,
    HarmonicVariable,
    _show_ansatz,
    _coordinate_transform,
    @eqtest,
    _create_harmonic_variable,
    get_variables_nums

# Setup
@variables t T
nat_var = declare_variable("x", t)
rot_var = declare_variable("u", T)
ω = declare_variable("ω")

@testset "Constructor" begin
    hv = HarmonicVariable(rot_var, "test_name", "u", ω, nat_var)
    @eqtest hv.symbol == rot_var
    @test hv.name == "test_name"
    @test hv.type == "u"
    @eqtest hv.ω == ω
    @eqtest hv.natural_variable == nat_var
end

@testset "Show and Display" begin
    hv = HarmonicVariable(rot_var, "test_name", "u", ω, nat_var)
    # Test show output
    io = IOBuffer()
    show(io, hv)
    @test String(take!(io)) == "Harmonic variable u(T) for harmonic ω of x(t)\n"

    # Test _show_ansatz
    @test contains(_show_ansatz(hv), "cos(ωt)")

    @variables x
    hv = HarmonicVariable(x)
    @test repr(hv) == "Harmonic variable x\n"

    @test _show_ansatz(hv) == "x"
end

@testset "Coordinate Transforms" begin
    new_var = declare_variable("y")
    # Test u-type transform
    @eqtest _coordinate_transform(new_var, ω, t, "u") == new_var * cos(ω * t)
    # Test v-type transform
    @eqtest _coordinate_transform(new_var, ω, t, "v") == new_var * sin(ω * t)
    # Test a-type transform
    @eqtest _coordinate_transform(new_var, ω, t, "a") == new_var
end

@testset "Create Harmonic Variable" begin
    rule, hvar = _create_harmonic_variable(nat_var, ω, t, "u"; new_symbol="a")
    a = declare_variable("a", t)
    @eqtest rule == a * cos(ω * t)
    @test hvar isa HarmonicVariable
    @test rule isa Num
    @test hvar.type == "u"
    @eqtest hvar.ω == ω
end

@testset "Substitution" begin
    hv = HarmonicVariable(nat_var, "test_name", "u", ω, nat_var)
    new_var = declare_variable("z")
    rules = Dict(hv => new_var)

    # Test equation substitution
    eq = nat_var^2 + ω
    result = substitute_all(eq, rules)
    @eqtest result == new_var^2 + ω

    # Test number substitution
    new_hv = substitute_all(hv, Dict(ω => 2))
    @test new_hv.ω == 2

    # Test variable substitution
    new_hv = substitute_all(hv, Dict(ω => new_var))
    @test isequal(new_hv.ω, new_var)
end
