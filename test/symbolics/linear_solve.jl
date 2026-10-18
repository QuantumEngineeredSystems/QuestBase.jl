using Test
using Symbolics
using LinearAlgebra: LinearAlgebra, det
using Random: Random
using QuestBase: QuestBase, @eqtest

@testset "fraction-free linear solve" begin
    using LinearAlgebra
    using QuestBase: fraction_free_linear_solve

    @variables a b c x y z
    equations = [2x + y - z ~ a, x + 3y + z ~ b, 2x - y + 4z ~ c]
    solution = fraction_free_linear_solve(equations, [x, y, z])
    @eqtest solution[1] == (13a - 3b + 4c) / 31
    @eqtest solution[2] == (-2a + 10b - 3c) / 31
    @eqtest solution[3] == (-7a + 4b + 5c) / 31

    # A zero diagonal requires a row pivot.
    pivoted = fraction_free_linear_solve([y ~ a, x + y ~ b], [x, y])
    @eqtest pivoted[1] == b - a
    @eqtest pivoted[2] == a

    @test_throws LinearAlgebra.SingularException fraction_free_linear_solve(
        [x + y ~ a, 2x + 2y ~ b], [x, y]
    )

    # Compare arbitrary nonsingular integer systems against exact rational arithmetic.
    rng = Random.MersenneTwister(0x5eed)
    for dimension in 2:5
        matrix = rand(rng, -4:4, dimension, dimension)
        while iszero(det(matrix))
            matrix = rand(rng, -4:4, dimension, dimension)
        end
        rhs = rand(rng, -4:4, dimension)
        symbolic_variables = only(@variables q[1:dimension])
        random_equations = [
            sum(
                matrix[row, column] * symbolic_variables[column] for column in 1:dimension
            ) ~ rhs[row] for row in 1:dimension
        ]
        actual = fraction_free_linear_solve(random_equations, symbolic_variables)
        expected = Rational{Int}.(matrix) \ Rational{Int}.(rhs)
        for index in eachindex(actual)
            @eqtest actual[index] == expected[index]
        end
    end
end

@testset "fraction-free linear solve splits independent blocks" begin
    using LinearAlgebra
    using QuestBase: fraction_free_linear_solve, system_blocks

    @variables a b c x y z w p q

    # Two 2x2 blocks, interleaved so the decomposition cannot rely on contiguous indices.
    # Each solution carries only its own block's determinant: solved whole, every one of
    # them would come out over the product of both.
    equations = [p * x + q * z ~ a, x - z ~ b, p * y + w ~ c, y - w ~ 0]
    solution = fraction_free_linear_solve(equations, [x, y, z, w])
    @eqtest solution[1] == (a + b * q) / (p + q)
    @eqtest solution[3] == (a - b * p) / (p + q)
    @eqtest solution[2] == c / (1 + p)
    @eqtest solution[4] == c / (1 + p)

    # Transcendental coefficients do not change the decomposition, and the trig-free block
    # is not dragged through the rotation block's determinant.
    trig = fraction_free_linear_solve(
        [cos(p) * x - sin(p) * z ~ a, sin(p) * x + cos(p) * z ~ b, y ~ c], [x, y, z]
    )
    @eqtest trig[2] == c
    # Put the rotation on a rational point of the unit circle so the residual is exact.
    numeric = Dict(cos(p) => 3 // 5, sin(p) => 4 // 5, a => 1, b => 2)
    @eqtest substitute(cos(p) * trig[1] - sin(p) * trig[3], numeric) == 1
    @eqtest substitute(sin(p) * trig[1] + cos(p) * trig[3], numeric) == 2

    # Whole matrix nonzero: one component, so there is nothing to split and the solve runs
    # as before.
    @test system_blocks(Num[1 1; 1 1]) == [([1, 2], [1, 2])]

    # A structurally singular pattern has an unbalanced component, and is left to the full
    # solve to reject rather than being silently split.
    @test system_blocks(Num[1 1; 0 0]) == [([1, 2], [1, 2])]
    @test_throws LinearAlgebra.SingularException fraction_free_linear_solve(
        [x + y ~ a, 0 * x + 0 * y ~ b], [x, y]
    )

    # Random block-diagonal integer systems against exact rational arithmetic.
    rng = Random.MersenneTwister(0xb10c)
    for _ in 1:8
        sizes = rand(rng, 1:3, 3)
        dimension = sum(sizes)
        matrix = zeros(Int, dimension, dimension)
        offset = 0
        for size in sizes
            block = rand(rng, -4:4, size, size)
            while iszero(det(block))
                block = rand(rng, -4:4, size, size)
            end
            matrix[(offset + 1):(offset + size), (offset + 1):(offset + size)] = block
            offset += size
        end
        permutation = Random.randperm(rng, dimension)
        matrix = matrix[permutation, :]
        rhs = rand(rng, -4:4, dimension)
        symbolic_variables = only(@variables r[1:dimension])
        random_equations = [
            sum(
                matrix[row, column] * symbolic_variables[column] for column in 1:dimension
            ) ~ rhs[row] for row in 1:dimension
        ]
        actual = fraction_free_linear_solve(random_equations, symbolic_variables)
        expected = Rational{Int}.(matrix) \ Rational{Int}.(rhs)
        for index in eachindex(actual)
            @eqtest actual[index] == expected[index]
        end
    end
end

@testset "rearrange! reduces determinant denominators" begin
    # Fraction-free elimination produces one determinant denominator instead of nested
    # LU fractions. For a trigonometric ansatz that determinant contains
    # cos² + sin² = 1, which `rearrange!` must collapse before downstream averaging.
    using QuestBase: HarmonicEquation, HarmonicVariable, DifferentialEquation
    using QuestBase: CosineQuadrature, SineQuadrature
    using QuestBase: rearrange!, d
    @variables ω t T u1(T) v1(T) x(t)

    eqs = [
        cos(ω * t) * d(u1, T) + sin(ω * t) * d(v1, T) ~ u1,
        -sin(ω * t) * d(u1, T) + cos(ω * t) * d(v1, T) ~ v1,
    ]
    hvars = [
        HarmonicVariable(u1, "u", CosineQuadrature, ω, x),
        HarmonicVariable(v1, "v", SineQuadrature, ω, x),
    ]
    natural_eq = DifferentialEquation(d(x, t, 2) + x ~ 0, x)
    eom = HarmonicEquation(eqs, hvars, natural_eq)

    rearrange!(eom, d([u1, v1], T))
    lhs = QuestBase.Num.(getfield.(eom.equations, :lhs))
    @eqtest lhs[1] == u1 * cos(ω * t) - v1 * sin(ω * t)
    @eqtest lhs[2] == u1 * sin(ω * t) + v1 * cos(ω * t)
end
