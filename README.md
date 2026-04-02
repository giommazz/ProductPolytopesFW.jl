# ProductPolytopesFW

Julia package and experiment scripts accompanying the AISTATS 2026 paper on the linear convergence of Frank-Wolfe algorithms on product polytope.

## Associated paper

This repository contains the code accompanying the paper

> Iommazzo, G., Martinez-Rubio, D., Criado, F., Wirth, E., and Pokutta, S.  
> **Linear Convergence of the Frank-Wolfe Algorithm over Product Polytopes.**  
> Accepted at *AISTATS 2026*.  
> Available as preprint at [arXiv:2505.11259](https://arxiv.org/abs/2505.11259).

The repository includes the implementation used for the computational results in the paper, together with additional experiment and benchmarking scripts.

## Citation

If you use this repository in your research, please cite:

```bibtex
@article{IMC+25,
  title={Linear Convergence of the Frank-Wolfe Algorithm over Product Polytopes},
  author={Iommazzo, Gabriele and Mart{\\'i}nez-Rubio, David and Criado, Francisco and Wirth, Elias and Pokutta, Sebastian},
  journal={arXiv preprint arXiv:2505.11259},
  year={2025},
  doi={10.48550/arXiv.2505.11259},
  url={https://arxiv.org/abs/2505.11259}
}
```

## Setup

Install as an unregistered Julia package:

```julia
using Pkg
Pkg.add(url="https://github.com/giommazz/ProductPolytopesFW.jl")
```

Then load it with:

```julia
using ProductPolytopesFW
```

For running experiments from the repository source:

```bash
cd /path/to/ProductPolytopesFW.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

The repository currently targets Julia `1.12` (`Project.toml` compat) and was verified locally with Julia `1.12.5`.

## Configuration

All experiment parameters live in `examples/config.yml`.
The example compute/plot scripts read this file directly; `examples/slurm_loop.jl` generates per-run config copies automatically.

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
- `convex_hull_backend: "matrix"`(this package's `ProductPolytopesFW.MatrixConvexHullLMO`) vs `"vector"` (`FrankWolfe.ConvexHullLMO`)
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
Run experiments on Slurm via

- `examples/slurm_experiments.sh`: submit one run with explicit script + results folder + config.
- `examples/slurm_loop.jl`: generate per-run configs/scripts and submit a `(k, n)` grid from CLI arguments.

`examples/slurm_experiments.sh` derives the repository root from its own location, and the cluster-specific `#SBATCH` lines for partition, node pinning, and email notifications are commented out by default. Enable only the directives that your cluster actually supports.

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
- `slurm_loop.jl` writes generated scripts/configs under `<results_dir>/slurm_generated/`.
- Set `JULIA_BIN=/path/to/julia` if the cluster does not expose `julia` in `PATH`.
- Optional `#SBATCH` directives in `examples/slurm_experiments.sh` are commented out by default; edit only what your cluster requires.

## Testing

Run the full deterministic test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Optional direct invocation:

```bash
julia --project=. test/runtests.jl
```
