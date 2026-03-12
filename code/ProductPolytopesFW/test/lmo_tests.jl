@testset "Matrix/vector LMO behavior" begin
    V = [
         1.0  0.0 -1.0;
        -2.0  1.0  0.5;
         0.3 -0.7  2.0;
         1.1  1.2  1.3;
         0.1  0.2  0.3;
    ]
    direction = [0.2, -1.5, 0.9]

    @testset "MatrixConvexHullLMO extreme point matches naive argmin" begin
        scores = vec(V * direction)
        idx = argmin(scores)
        expected = collect(@view V[idx, :])

        lmo = ProductPolytopesFW.MatrixConvexHullLMO(V; cache_cap=2, use_optimized_search=true)
        s = FrankWolfe.compute_extreme_point(lmo, direction)

        @test s isa Vector{Float64}
        @test s == expected
    end

    @testset "Optimized search and row-scan return same extreme point" begin
        lmo_opt = ProductPolytopesFW.MatrixConvexHullLMO(V; cache_cap=2, use_optimized_search=true)
        lmo_scan = ProductPolytopesFW.MatrixConvexHullLMO(V; cache_cap=2, use_optimized_search=false)

        s_opt = FrankWolfe.compute_extreme_point(lmo_opt, direction)
        s_scan = FrankWolfe.compute_extreme_point(lmo_scan, direction)

        @test s_opt == s_scan
    end

    @testset "Cache on/off parity and dense output" begin
        lmo_cache = ProductPolytopesFW.MatrixConvexHullLMO(V; cache_cap=3, use_optimized_search=true)
        lmo_nocache = ProductPolytopesFW.MatrixConvexHullLMO(V; cache_cap=0, use_optimized_search=true)

        dirs = (
            [0.2, -1.5, 0.9],
            [-0.1, 0.3, 2.0],
            [1.0, 1.0, 1.0],
            [0.2, -1.5, 0.9],
        )
        for d in dirs
            s_cache = FrankWolfe.compute_extreme_point(lmo_cache, d)
            s_nocache = FrankWolfe.compute_extreme_point(lmo_nocache, d)

            @test s_cache == s_nocache
            @test s_cache isa Vector{Float64}

            out = zeros(3)
            returned = FrankWolfe.compute_extreme_point(lmo_cache, d; v=out)
            @test returned === out
            @test returned isa Vector{Float64}
            @test returned == s_cache
        end
    end
end
