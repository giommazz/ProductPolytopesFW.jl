# ProductPolytopesAFW

Julia package and experiment scripts for Frank-Wolfe variants on product polytope feasibility instances.

## Setup

```bash
cd code/ProductPolytopesAFW
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

The repository currently targets Julia `1.12` (`Project.toml` compat).

## Configuration

All experiment parameters live in `examples/config.yml`.

## Current experiment scripts

### Point-cloud workflow

Runs experiments on polytopes generated from sampled point clouds and solved with convex-hull LMOs.
Supports multiple instance-generation and intersection settings via `examples/config.yml`: modify this file to generate+solve different instances

```bash
mkdir -p examples/results_linesearch_point_clouds/terminal_logs
julia --project=. examples/compute_intersection_point_clouds_full_warmup.jl 2>&1 | tee examples/results_linesearch_point_clouds/terminal_logs/run_$(date +%Y%m%d_%H%M%S).log
```

### Vertex-facet workflow

Runs experiments on polytope instances with vertex-facet geometry, where one vertex of a generalized `\ell_1` ball (diamond) touches a facet of a generalized `\ell_\infty` box.
Modify main parameters in `compute_intersection_vertex_facet_full_warmup.jl` to generate+solve different instances

```bash
mkdir -p examples/results_linesearch_vertex_facet/terminal_logs
julia --project=. examples/compute_intersection_vertex_facet_full_warmup.jl 2>&1 | tee examples/results_linesearch_vertex_facet/terminal_logs/run_$(date +%Y%m%d_%H%M%S).log
```

### Plotting from saved logs

```bash
julia --project=. examples/plot_from_logs.jl
julia --project=. examples/plot_from_logs_avg.jl
```

## Profiling and backend benchmarking

To profile wall time and memory of one run:

```bash
/usr/bin/time -v julia --project=. examples/compute_intersection_point_clouds_full_warmup.jl 2>&1 | tee examples/results_bench_cvxh_backend/terminal_example.log
```

To compare convex-hull backends, run the same script while changing in `examples/config.yml`:
- `convex_hull_backend: "matrix"`(this package's `ProductPolytopesAFW.MatrixConvexHullLMO`) vs `"vector"` (`FrankWolfe.ConvexHullLMO`)
- `matrix_lmo_use_optimized_search: true/false` (only for matrix backend)
- `matrix_lmo_cache_cap: -1/0/>0` (only for matrix backend)

Then extract a compact summary from terminal logs:

```bash
for f in examples/results_bench_cvxh_backend/terminal_*.log; do
  echo "== $f =="
  grep -E "Elapsed time main:|Maximum resident set size|Exit status:" "$f"
done
```

## SLURM
Run experiments on slurm via

- `examples/slurm_experiments.sh`: submit one run with explicit script + results folder + config.
- `examples/slurm_loop.jl`: generate per-run configs/scripts and submit a `(k, n)` grid.

Single submission for point-cloud-type instances:

```bash
sbatch examples/slurm_experiments.sh examples/compute_intersection_point_clouds_full_warmup.jl examples/results_linesearch_point_clouds examples/config.yml
```

Single submission for vertex-facet-type instances:

```bash
sbatch examples/slurm_experiments.sh examples/compute_intersection_vertex_facet_full_warmup.jl examples/results_linesearch_vertex_facet examples/config.yml
```

Grid submission example command for point-cloud-type instances (with dry-run first):

```bash
julia --project=. examples/slurm_loop.jl examples/compute_intersection_point_clouds_full_warmup.jl examples/results_linesearch_point_clouds examples/config.yml 2,3 103,207 555 --dry-run
julia --project=. examples/slurm_loop.jl examples/compute_intersection_point_clouds_full_warmup.jl examples/results_linesearch_point_clouds examples/config.yml 2,3 103,207 555
```

Notes:
- `slurm_loop.jl` overrides `--cpus-per-task` per job as `k`.
- Set `JULIA_BIN=/path/to/julia` if the cluster does not expose `julia` in `PATH`.
- Adjust `SBATCH` directives in `examples/slurm_experiments.sh` to match your cluster.

## Testing

Run the full deterministic test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Optional direct invocation:

```bash
julia --project=. test/runtests.jl
```
