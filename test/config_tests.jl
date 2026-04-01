@testset "Config parsing and updates" begin
    @testset "Config YAML parsing includes new fields" begin
        yaml_text = """
        k: 2
        n: 5
        n_points: [7, 8]
        ub_n_points: 50
        vertex_sampling: "ellipsoid"
        sphere_radius_factor: 0.8
        target_tolerance: 1.0e-7
        target_tolerance_opt: 1.0e-9
        max_iterations: 120
        max_iterations_opt: 240
        max_print_iterations: 20
        seed: 123
        verbose: false
        cvxhflag: true
        convex_hull_backend: "vector"
        matrix_lmo_cache_cap: 0
        matrix_lmo_use_optimized_search: false
        intersection:
          anchor: "p1_p2_segment"
          anchor_t: 0.25
          reference_point: "vertex"
        stepsize_strategy: 1
        """

        cfg_path = tempname() * ".yml"
        open(cfg_path, "w") do io
            write(io, yaml_text)
        end

        cfg = ProductPolytopesFW.Config(cfg_path)
        @test cfg.convex_hull_backend == "vector"
        @test cfg.matrix_lmo_cache_cap == 0
        @test cfg.matrix_lmo_use_optimized_search == false
        @test cfg.intersection_anchor == "p1_p2_segment"
        @test cfg.intersection_reference_point == "vertex"
        @test cfg.intersection_anchor_t == 0.25
        @test cfg.n_points == [7, 8]
    end

    @testset "modify_config updates requested fields and preserves others" begin
        base = ProductPolytopesFW.Config()
        updated = ProductPolytopesFW.modify_config(
            base;
            max_iterations=222,
            max_print_iterations=33,
            convex_hull_backend="vector",
            matrix_lmo_cache_cap=12,
            matrix_lmo_use_optimized_search=false,
            intersection_anchor="global_mean",
            intersection_reference_point="vertex",
            intersection_anchor_t=0.7,
        )

        @test updated.max_iterations == 222
        @test updated.max_print_iterations == 33
        @test updated.convex_hull_backend == "vector"
        @test updated.matrix_lmo_cache_cap == 12
        @test updated.matrix_lmo_use_optimized_search == false
        @test updated.intersection_anchor == "global_mean"
        @test updated.intersection_reference_point == "vertex"
        @test updated.intersection_anchor_t == 0.7
        @test updated.k == base.k
        @test updated.n == base.n
        @test updated.stepsize_strategy == base.stepsize_strategy
    end
end
