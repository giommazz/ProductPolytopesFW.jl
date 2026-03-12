@testset "Objective and gradient parity" begin
    @testset "convex_feasibility_objective_v1 ≈ convex_feasibility_objective_v2b" begin
        rng = MersenneTwister(20260306)
        x = ProductPolytopesFW.random_blockvector(3, 7; rng=rng)
        f1 = ProductPolytopesFW.convex_feasibility_objective_v1(x)
        f2 = ProductPolytopesFW.convex_feasibility_objective_v2b(x)
        @test isapprox(f1, f2; rtol=1e-12, atol=1e-12)
    end

    @testset "convex_feasibility_gradient_v1! ≈ convex_feasibility_gradient_v2!" begin
        rng = MersenneTwister(314159)
        x = ProductPolytopesFW.random_blockvector(4, 5; rng=rng)
        g1 = ProductPolytopesFW.zeros_blockvector(4, 5)
        g2 = ProductPolytopesFW.zeros_blockvector(4, 5)

        ProductPolytopesFW.convex_feasibility_gradient_v1!(g1, x)
        ProductPolytopesFW.convex_feasibility_gradient_v2!(g2, x)

        for i in eachindex(g1.blocks)
            @test isapprox(g1.blocks[i], g2.blocks[i]; rtol=1e-12, atol=1e-12)
        end
    end
end
