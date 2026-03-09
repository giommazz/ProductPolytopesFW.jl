# Benchmarking `ConvexHullLMO` against `MatrixConvexHullLMO`

Compare the following variants:

- `vec_best`: `ConvexHullLMO` with most efficient vector-based setup
- `mat_opt_cacheauto`: `MatrixConvexHullLMO` with optimized matrix-vector scan and default cache
- `mat_opt_cacheoff`: same optimized scan, cache disabled
- `mat_scan_cacheauto`: `MatrixConvexHullLMO` with row-by-row scan and default cache
- `mat_scan_cacheoff`: same row-by-row scan, cache disabled

## Tables
Summary tables are as follows.

### Oracle table (`results/oracle_table`)
Measures only repeated calls to `compute_extreme_point`, without running Frank-Wolfe.

For each case, the total number of oracle calls is

`oracle_calls = oracle_repetitions * k * direction_batch_size`

where:

- `k` is the number of blocks / LMOs in the product
- `direction_batch_size` is the number of random directions tested per block in each repetition
- a batch means one full pass over all `k` blocks, using `direction_batch_size` random directions for each block
- `oracle_repetitions` is the number of times this same batch structure is repeated


`oracle_table` uses the following config: `oracle_repetitions = 200`, `k = 2`, and `direction_batch_size = 32` so each oracle case performs `200 * 2 * 32 = 12800` calls to `compute_extreme_point`.

### Algorithm table (`results/algorithm_table`)
Runs a very short solver execution with `FrankWolfe.frank_wolfe` and `FrankWolfe.away_frank_wolfe`.

Both use `SafeGoldenratio` as line search. Goal: measure backend cost during algorithm execution.

### Meaning of the columns

Both tables use the same main columns:

- `case_id`: benchmark variant name
- `n`: ambient dimension
- `n_points`: number of atoms / vertices in each point cloud
- `elapsed_s`: total wall-clock time of the whole workload
- `alloc_bytes`: total number of bytes allocated by Julia during the whole workload
- `alloc_count`: total number of Julia allocation events during the whole workload
- `max_rss_kb`: peak memory of the whole process, in KB

The call-count column is different:

- `oracle_calls` in `oracle_table`: total number of `compute_extreme_point` calls, see above
- `iterations` in `algorithm_table`: total number of FW or AFW iterations

Instead, `n_points = n_points_multiplier * n` and, with the default `n_points_multiplier = 2` one has:

- `n = 100  -> n_points = 200`
- `n = 500  -> n_points = 1000`
- `n = 1000 -> n_points = 2000`
- `n = 2000 -> n_points = 4000`

### Meaning of the five configurations

- `vec_best`: baseline `ConvexHullLMO`, built from a *vector* of row-atoms
- `mat_opt_*`: matrix backend using one *dense matrix-vector product* (`mul!` / BLAS-backed) plus an `argmin` over the resulting scores
- `mat_scan_*`: matrix backend using a row-by-row scan (one dot product per row)
- `cacheauto`: materialized-vertex cache enabled with the backend default policy
- `cacheoff`: cache disabled (`cache_cap = 0`)

`vec_best` VS `mat_scan` both score one candidate atom at a time and keep the best one found so far. The difference is in how the atoms are stored and accessed:

- `vec_best` uses the standard `ConvexHullLMO` on a Julia collection of row-atoms.
- `mat_scan` uses the matrix-backed LMO, but still evaluates the rows one by one: it is the same style of search (one score at a time), but with the atoms stored in one dense matrix instead of a generic vector-like container.

`mat_opt` changes the implementation: it computes all row scores at once through one dense matrix-vector multiplication and only then takes the min score (exploits contiguous storage and BLAS routines).

The cache is used only by the matrix backend. It stores dense copies of vertices that were selected recently by the oracle. This matters because the matrix backend stores atoms as matrix rows. When one row is selected, the backend often needs to materialize that row as a standalone dense `Vector`. If later the same row index is selected again, `cacheauto` lets the code reuse the already materialized vector instead of copying that row from the matrix again.

So the cache is not a cache of the whole active set and not a cache of arbitrary directions. It is a cache of recently returned atoms, keyed by vertex row index. In practice, useful when the same atoms are selected repeatedly across nearby oracle calls.

`cacheauto` means the cache size is chosen automatically by the matrix LMO using a heuristic rule. `cacheoff` means the cache is disabled completely, so every selected row must be materialized again when needed.

## How time and memory were captured

Inside Julia, the scripts measure:

- time with `@elapsed`
- allocated bytes with `@allocated`
- allocation count with `@allocations`

The measurement is taken around the whole benchmark workload, with `GC.gc()` calls before each macro to reduce noise.

At the process level, `run_suite.sh` wraps each Julia command with:

- `/usr/bin/time -v`

and extracts:

- `Maximum resident set size (kbytes)`

from the raw logs. This is what appears as `max_rss_kb` in the tables.

Interpretation:

- `alloc_bytes`:cumulative allocation traffic over the whole run
- `alloc_count`: cumulative number of allocations over the whole run
- `max_rss_kb`: peak RAM footprint at one instant

## How to regenerate the tables

Run the full suite:

`bash benchmarks_local_cvxh_lmo/run_suite.sh`

or regenerate only the tables from existing raw logs and per-case CSV files:

`julia --project=. benchmarks_local_cvxh_lmo/summarize_results.jl benchmarks_local_cvxh_lmo/config.yml`
