using Test
using Symbolics
using LinearAlgebra

using QuestBase: strip_derivative, dummy_symbolic_Jacobian, is_identity, hasnan, d, @eqtest

@testset "strip_derivative" begin
    @variables t x(t) y(t)

    @eqtest strip_derivative(x) == x
    @eqtest strip_derivative(d(x, t)) == x
    @eqtest strip_derivative(d(d(x, t), t)) == x

    # should leave non-derivative expressions untouched
    @eqtest strip_derivative(x + y) == x + y
end

@testset "eqtest_equal containers and complex" begin
    @variables a b c

    @eqtest [a, b] == [a, b]
    @test !QuestBase._eqtest_equal([a], [a, b])

    @eqtest (a, b, c) == (a, b, c)
    @test !QuestBase._eqtest_equal((a, b), (a, c))

    z = a + im * b
    @eqtest z == z
    @eqtest (a + 0im) == a
    @eqtest a == (a + 0im)
end

@testset "identity and NaN helpers" begin
    @variables a

    I2 = Matrix{Num}(LinearAlgebra.I, 2, 2)
    @test is_identity(I2)

    A = Num[1 0; 0 (1 + a)]
    @test !is_identity(A)

    J = dummy_symbolic_Jacobian(3)
    @test hasnan(J)
    @test !hasnan(I2)
    @test !is_identity(J)
end
