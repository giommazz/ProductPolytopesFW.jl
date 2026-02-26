# TODOS

Next tasks to implement or explore to generate intersection instances where AFW linear behavior is visible and where FW vs AFW differences are easier to interpret.

## Edge-Biased Sampling for `box_uniform` (`box_boundary_power`)
### Why
`ConvexHullLMO` scans all sampled points, but only extreme points can be selected by the oracle. With uniform interior sampling, many points are redundant and increase runtime.
### Idea
Keep `vertex_sampling: "box_uniform"` and add a parameter that biases coordinates toward LB/UB. This should increase the fraction of useful boundary points.
### Implementation sketch
1. Add `box_boundary_power: p` to `examples/config.yml` with `p >= 1.0` (`1.0` keeps current behavior).
2. In `generate_polytope()` for `vertex_sampling == "box_uniform"`:
   - sample `u ~ Uniform(0,1)`,
   - with probability `0.5` set `t = u^p`, otherwise `t = 1 - u^p`,
   - map each coordinate as `x = LB + (UB - LB) * t`.
3. Optional: add `box_boundary_eps` and clamp `t` to `[eps, 1 - eps]` to keep points strictly inside the box.
4. Update `generate_filename` to include the bias parameter (example: `box-unif-p2.0`).
5. Potential files: `examples/config.yml`, `src/config.jl`, `src/polytope_generation.jl`, `src/utils.jl`.



## Regularize Box Shape (`box_width_min_factor`)
### Why
For `sphere`, radius depends on `min(h[d])`, so a single narrow coordinate can collapse the sampled geometry. For `ellipsoid`, axes directly inherit box anisotropy and can become highly squashed.
### Idea
Control coordinate-wise box widths by enforcing a minimum width factor to avoid very small `min(h)` and extreme `max/min` ratios.
### Implementation sketch
1. Add `box_width_min_factor: wmin` to `examples/config.yml` with `0 < wmin <= 1`.
2. In `generate_nonintersecting_bounds()`, enforce `(UB - LB) >= wmin * stepsize` per coordinate (for example, sample widths in `[wmin * stepsize, stepsize]`).
3. Optional: keep/extend diagnostics with half-length ratio statistics.
4. Update filename metadata when `wmin` is active.
5. Potential files: `examples/config.yml`, `src/config.jl`, `src/polytope_generation.jl`, `src/utils.jl`.



## Non-Monotone Box Separation (Harder Non-Intersecting for `sphere`/`ellipsoid`)
### Why
Current generation separates boxes in all coordinates. Combined with small sampled polytopes, this can create trivial non-intersecting instances.
### Idea
Generate disjoint boxes without forcing coordinate-wise monotone separation, while preserving guaranteed disjointness.
### Implementation sketch
1. Add `bounds_separation: "monotone" | "random"` to `examples/config.yml`.
2. For `"random"`, sample box centers and enforce disjointness via a geometric rule (for example `L_inf` separation plus margin).
3. Add a cheap consistency check: disjoint boxes imply disjoint sampled polytopes (since points remain inside boxes).
4. Tag filenames with the separation mode (`sep=mono` or `sep=rand`).
5. Potential files: `examples/config.yml`, `src/config.jl`, `src/polytope_generation.jl`, `src/utils.jl`.



## Global Preconditioning (`none`/`diag`/`whiten`)
### Why
Anisotropy and heterogeneous coordinate scales can hurt conditioning and slow FW variants. A shared global transformation can improve geometry.
### Idea
Apply the same affine transform to all polytopes:
1. `diag`: subtract global mean and scale each coordinate by standard deviation (or range)
2. `whiten`: apply whitening with regularization.
### Implementation sketch
1. Add `preconditioning: "none" | "diag" | "whiten"` and supporting parameters (for example `whiten_eps`) to `examples/config.yml`.
2. Implement helper functions in `src/polytope_utils.jl` to compute transformation parameters.
3. In `generate_polytopes()`, transform all generated/shifted vertices with the same map.
4. If structured LMOs are used later, handle compatibility explicitly (`diag` is straightforward on box bounds; `whiten` may require restrictions or extra handling).
5. Add preconditioning tags to filenames (`pc=diag`, `pc=whiten`).
6. Potential files: `examples/config.yml`, `src/config.jl`, `src/polytope_utils.jl`, `src/polytope_generation.jl`, `src/utils.jl`.



## Deeper `p1_random` Anchor (Dirichlet Alpha > 1)
### Why
Current `p1_random` anchor may place too much weight on a few vertices, producing anchors close to faces and worsening local geometry.
### Idea
Use Dirichlet-distributed convex weights with `alpha > 1` to produce more balanced combinations and anchors deeper in `P1`.
### Implementation sketch
1. Add `p1_random_dirichlet_alpha: 1.0` to `examples/config.yml` (`1.0` preserves current behavior).
2. Implement Dirichlet sampling via gamma draws (`w_i ~ Gamma(alpha, 1)` then normalize).
3. Use these weights in `random_convex_combination()`.
4. Add `diralpha` to filenames/metadata.
5. Potential files: `examples/config.yml`, `src/config.jl`, `src/polytope_utils.jl`, `src/utils.jl`.



## Vertex-Facet Instances (Sebastian Suggestion) with Structured LMOs
### Why
FW often identifies the active/optimal face asymptotically, while AFW can identify it in finite time under suitable
### Idea
Design instances to expose support/face identification differences. For example, for `k = 2`:
- `P1` is a box/hypercube
- `P2` is a simplex (small convex hull) with one vertex touching the relative interior of facet of `P1`, and other vertices outside along the facet normal.
### Implementation sketch
1. Extend config to declare polytope type per block (example: `polytope_types: ["box", "cvxhull"]`).
2. Add a dedicated generator family for this construction.
3. Integrate it with current experiment scripts (`compute_intersection_custom_full_warmup`, SLURM scripts).
4. Distinguish this family in plot/file naming.
5. Potential files: `examples/config.yml`, `src/config.jl`, `src/lmo_utils.jl`, `src/polytope_generation.jl`, `src/utils.jl`, `examples/compute_intersection_custom_full_warmup.jl`.



## Support Mixed Legacy/Convex-Hull Polytopes
### Why
Legacy LMOs can reduce runtime wrt ConvexHullLMO
### Idea
Allow each block/polytope to be generated either as a sampled convex hull or as a structured object (box, simplex, L1-ball, and similar).
### Implementation sketch
1. Extend `Config` with per-block polytope type and parameters.
2. Update LMO creation logic to instantiate the correct oracle for each block.
3. Skip point-cloud generation for structured blocks.
4. Potential files: `examples/config.yml`, `src/config.jl`, `src/lmo_utils.jl`, `src/polytope_generation.jl`, `src/ProductPolytopesAFW.jl`.



## Vertex Pruning for Faster `ConvexHullLMO`
### Why
Even with better sampling, `n_points` can remain very large. Pruning can keep runtime manageable while preserving most useful extreme points.
### Idea
Select extreme points along random directions and keep only their union as a reduced candidate set.
### Implementation sketch
1. Add `vertex_pruning: "none" | "random_directions"` and `num_directions: M` to `examples/config.yml`.
2. Implement pruning in `src/polytope_utils.jl` by collecting `argmin`/`argmax` points for sampled directions.
3. Apply pruning after generation and before oracle construction.
4. Tag filenames with pruning settings (`prune=M`).
5. Potential files: `examples/config.yml`, `src/config.jl`, `src/polytope_utils.jl`, `src/polytope_generation.jl`, `src/utils.jl`.



## Improve Filenames and Metadata for Reproducibility
### Why
Many parameters materially affect results, check that all contribute to metadata
### Idea
Include compact but informative run metadata in filenames and saved artifacts.
### Implementation sketch
1. Encode effective `n_points` summary in filename (for example `npmin/npmax/npmean`).
2. Save full resolved config and geometry diagnostics in `.jld2`.
3. Read metadata in plotting scripts when available.
4. Potential files: `src/utils.jl`, `src/polytope_generation.jl`, `src/plotting_utils.jl`, `examples/plot_from_logs.jl`.



## Make Geometry Diagnostics Optional and Uniform
### Why
Current debug prints are useful locally but noisy in SLURM runs. A controlled and uniform diagnostics path is needed.
### Idea
Add a config flag and centralize diagnostic computation/printing.
### Implementation sketch
1. Add `print_geometry_stats: true | false` to `examples/config.yml`.
2. Move/centralize diagnostic helpers in `src/polytope_utils.jl`.
3. Keep a consistent print format for `box_uniform`, `sphere`, and `ellipsoid`.
4. Potential files: `examples/config.yml`, `src/config.jl`, `src/polytope_utils.jl`, `src/polytope_generation.jl`.



## Fix `log_times` Aggregation of Cumulative Timestamps
### Why
If trajectory time entries are cumulative, summing them over iterations overestimates total runtime.
### Idea
Use the last timestamp as total time when data is cumulative; otherwise sum deltas.
### Implementation sketch
1. Confirm semantics of the trajectory time field in `product_algorithms.jl`.
2. Update `log_times()` in `src/utils.jl` to handle cumulative vs delta modes correctly.
3. Keep CSV output consistent with the selected interpretation.
4. Potential files: `src/utils.jl`, `examples/compute_intersection_custom_full_warmup.jl`.



## Clarify Plotted Gaps (Primal vs FW, Global vs Local)
### Why
Observed cases with FW gap below primal gap suggest potential mismatch in logged metrics (global vs block-local gap) or plotting logic.
### Idea
Audit logged quantities and make plotted definitions explicit and consistent across methods.
### Implementation sketch
1. Verify which gap is returned/logged for each solver variant.
2. Optionally log both local and global gap for diagnostics.
3. Ensure primal plots use `f(x_k) - f*` with `f*` computed from the designated optimal run.
4. Potential files: `src/product_algorithms.jl`, `src/plotting_utils.jl`, `src/utils.jl`, `examples/compute_intersection_custom_full_warmup.jl`.