@testset "Tiny FW/AFW smoke runs" begin
    base_cfg = ProductPolytopesAFW.Config()
    cfg = ProductPolytopesAFW.modify_config(
        base_cfg;
        k=2,
        n=3,
        n_points=[4, 4],
        target_tolerance=1e-9,
        target_tolerance_opt=1e-12,
        max_iterations=12,
        max_print_iterations=1000,
        cvxhflag=true,
        convex_hull_backend="matrix",
        matrix_lmo_cache_cap=2,
        matrix_lmo_use_optimized_search=true,
        stepsize_strategy=0,
        verbose=false,
    )

    V1 = [
        0.0 0.0 0.0;
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0;
    ]
    V2 = [
        2.0 2.0 2.0;
        3.0 2.0 2.0;
        2.0 3.0 2.0;
        2.0 2.0 3.0;
    ]

    lmos = ProductPolytopesAFW.create_lmos(cfg, [V1, V2])
    prod_lmo = ProductPolytopesAFW.create_product_lmo(lmos)

    x_fw, v_fw, primal_fw, gap_fw, traj_fw = ProductPolytopesAFW.run_FullFW(cfg, FrankWolfe.frank_wolfe, prod_lmo)
    @test isfinite(primal_fw)
    @test isfinite(gap_fw)
    @test !isempty(traj_fw)
    @test x_fw !== nothing
    @test v_fw !== nothing

    x_afw, v_afw, primal_afw, gap_afw, traj_afw = ProductPolytopesAFW.run_FullAFW(cfg, prod_lmo)
    @test isfinite(primal_afw)
    @test isfinite(gap_afw)
    @test !isempty(traj_afw)
    @test x_afw !== nothing
    @test v_afw !== nothing
end
